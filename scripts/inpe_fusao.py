#!/usr/bin/env python3
"""Imagens do INPE (CBERS / Resourcesat) para um imóvel: busca, corrige o georreferenciamento,
faz a fusão pancromática + multiespectral e gera GeoTIFFs, mapa da fazenda, pranchas de zoom
por alerta e ficha técnica, tudo num .zip.

Entrada (variáveis de ambiente, enviadas pelo app no modo técnico):
  PEDIDO  identificador do pedido (nome do .zip)
  ANOS    ex.: "2008, 2013, atual"
  DADOS   JSON {car, nome, mun, imovel: geometria GeoJSON, alertas: [{id, fonte, ano, ha, geom}]}
Saída: saida/<PEDIDO>.zip
"""
import datetime as dt
import io
import json
import math
import os
import re
import sys
import urllib.request
import zipfile

import cv2
import numpy as np
import rasterio
from affine import Affine
from PIL import Image, ImageDraw, ImageFont
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds, transform_geom
from rasterio.windows import from_bounds
from shapely.geometry import box, mapping, shape
import warnings
warnings.filterwarnings("ignore")
from shapely.ops import transform as shp_transform
from pyproj import Transformer

os.environ.update(GDAL_HTTP_MAX_RETRY="6", GDAL_HTTP_RETRY_DELAY="4", GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                  CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif", GDAL_HTTP_TIMEOUT="180", VSI_CACHE="TRUE",
                  GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", GDAL_HTTP_MULTIPLEX="YES", GDAL_HTTP_VERSION="2")
STAC = "https://data.inpe.br/bdc/stac/v1/search"
MARCO = dt.date(2008, 7, 22)


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- produtos INPE
# cada produto: coleção MS, bandas (verde, vermelho, NIR), resolução; pan opcional; L4 = ortorretificado
PRODUTOS = [
    dict(nome="CBERS-4A WPM", ms="CB4A-WPM-L4-DN-1", b=("BAND2", "BAND3", "BAND4"), res_ms=8, pan=("CB4A-WPM-L4-DN-1", "BAND0"), res_pan=2, l4=True,
         ini=dt.date(2019, 12, 29), fim=None, sensores="WPM pancromática 2 m + multiespectral 8 m"),
    dict(nome="CBERS-4 PAN", ms="CB4-PAN10M-L4-DN-1", b=("BAND2", "BAND3", "BAND4"), res_ms=10, pan=("CB4-PAN5M-L4-DN-1", "BAND1"), res_pan=5, l4=True,
         ini=dt.date(2014, 12, 9), fim=None, sensores="PAN pancromática 5 m + multiespectral 10 m"),
    dict(nome="CBERS-4 MUX", ms="CB4-MUX-L4-DN-1", b=("BAND6", "BAND7", "BAND8"), res_ms=20, pan=("CB4-PAN5M-L4-DN-1", "BAND1"), res_pan=5, l4=True,
         ini=dt.date(2014, 12, 9), fim=None, sensores="PAN pancromática 5 m + MUX multiespectral 20 m"),
    dict(nome="CBERS-2B HRC/CCD", ms="CB2B-CCD-L2-DN-1", b=("CCD1XS_BAND2", "CCD1XS_BAND3", "CCD1XS_BAND4"), res_ms=20, pan=("CB2B-HRC-L2-DN-1", "BAND1"), res_pan=2.5, l4=False,
         ini=dt.date(2007, 9, 25), fim=dt.date(2010, 3, 11), sensores="HRC pancromática 2,5 m + CCD multiespectral 20 m"),
    dict(nome="Resourcesat-1 LISS-III", ms="RESOURCESAT1_LISS3_L2_DN-1", b=("BAND2", "BAND3", "BAND4"), res_ms=23.5, pan=None, l4=False,
         ini=dt.date(2009, 9, 10), fim=dt.date(2013, 9, 12), sensores="LISS-III multiespectral 23,5 m (sem pancromática)"),
    dict(nome="CBERS-2B CCD", ms="CB2B-CCD-L2-DN-1", b=("CCD1XS_BAND2", "CCD1XS_BAND3", "CCD1XS_BAND4"), res_ms=20, pan=None, l4=False,
         ini=dt.date(2007, 9, 25), fim=dt.date(2010, 3, 11), sensores="CCD multiespectral 20 m (sem pancromática no dia)"),
    dict(nome="CBERS-2 CCD", ms="CB2-CCD-L2-DN-1", b=("CCD1XS_BAND2", "CCD1XS_BAND3", "CCD1XS_BAND4"), res_ms=20, pan=None, l4=False,
         ini=dt.date(2003, 10, 28), fim=dt.date(2009, 1, 7), sensores="CCD multiespectral 20 m"),
]
FAIXAS = {"verde": "0,52–0,59 µm", "vermelho": "0,63–0,69 µm", "NIR": "0,77–0,89 µm"}


def stac(colecao, bbox, ini, fim):
    corpo = json.dumps({"collections": [colecao], "bbox": bbox, "datetime": f"{ini}T00:00:00Z/{fim}T23:59:59Z", "limit": 200}).encode()
    for t in range(4):
        try:
            req = urllib.request.Request(STAC, corpo, {"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=120)).get("features", [])
        except Exception as e:
            log("   STAC falhou:", e)
    return []


def por_data(itens):
    g = {}
    for it in itens:
        g.setdefault(it["properties"]["datetime"][:10], []).append(it)
    return g


# ---------------------------------------------------------------- leitura
def ler(href, bb_wgs, pad):
    """Lê uma janela (bbox WGS84 + pad em graus) de um COG remoto. Devolve (array, transform, crs) ou None."""
    try:
        with rasterio.open("/vsicurl/" + href) as s:
            b = transform_bounds("EPSG:4326", s.crs, bb_wgs[0] - pad, bb_wgs[1] - pad, bb_wgs[2] + pad, bb_wgs[3] + pad)
            w = from_bounds(*b, s.transform).round_offsets().round_lengths()
            a = s.read(1, window=w, boundless=True, fill_value=0)
            if not (a > 0).any():
                return None
            return a, s.window_transform(w), s.crs
    except Exception as e:
        log("   leitura falhou:", href.split("/")[-1], e)
        return None


def ler_rapido(href, bb_wgs, px=200):
    """Leitura reduzida (usa as visões gerais do COG) só para avaliar nuvem."""
    try:
        with rasterio.open("/vsicurl/" + href) as s:
            b = transform_bounds("EPSG:4326", s.crs, *bb_wgs)
            w = from_bounds(*b, s.transform)
            return s.read(1, window=w, boundless=True, fill_value=0, out_shape=(px, px), resampling=Resampling.average)
    except Exception:
        return None


def deslocar(peca, dx, dy):
    a, t, c = peca
    return a, Affine(t.a, t.b, t.c - dx, t.d, t.e, t.f - dy), c


def mosaico(pecas, grade, rs=Resampling.cubic):
    crs, tr, W, H = grade
    out = np.zeros((H, W), np.float32)
    for a, t, c in pecas:
        tmp = np.zeros((H, W), np.float32)
        reproject(a.astype(np.float32), tmp, src_transform=t, src_crs=c, dst_transform=tr, dst_crs=crs,
                  resampling=rs, src_nodata=0, dst_nodata=0)
        m = (out == 0) & (tmp > 0)
        out[m] = tmp[m]
    return out


# ---------------------------------------------------------------- nuvem / escolha de cena
def nota_nuvem(item, prod, bb):
    """Nuvem e falta de dado sobre o imóvel, medidas na miniatura da cena (rápido: não baixa a imagem)."""
    try:
        th = item["assets"].get("thumbnail")
        if not th:
            return 50
        with rasterio.open("/vsicurl/" + item["assets"][prod["b"][1]]["href"]) as s:
            W, H, tr, crs = s.width, s.height, s.transform, s.crs
        req = urllib.request.Request(th["href"], headers={"User-Agent": "RA-AgroTec"})
        im = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(req, timeout=60).read())).convert("RGB")).astype(np.float32)
        b = transform_bounds("EPSG:4326", crs, *bb)
        inv = ~tr
        c0, r0 = inv * (b[0], b[3])
        c1, r1 = inv * (b[2], b[1])
        kx, ky = im.shape[1] / W, im.shape[0] / H
        x0, x1 = int(math.floor(min(c0, c1) * kx)) - 1, int(math.ceil(max(c0, c1) * kx)) + 1
        y0, y1 = int(math.floor(min(r0, r1) * ky)) - 1, int(math.ceil(max(r0, r1) * ky)) + 1
        if x1 <= 0 or y1 <= 0 or x0 >= im.shape[1] or y0 >= im.shape[0]:
            return 999
        rec = im[max(y0, 0):y1, max(x0, 0):x1]
        fora = 1 - rec.shape[0] * rec.shape[1] / max(1, (y1 - y0) * (x1 - x0))
        lum = rec.mean(axis=2)
        vazio = (lum < 3).mean()
        nuvem = (lum > 205).mean()
        return nuvem * 100 + (vazio + fora) * 100
    except Exception as e:
        log("   miniatura falhou:", e)
        return 50


def escolher(ano, bb, bb_folga):
    """Devolve (produto, data, itens_ms, itens_pan, nota) da melhor cena para o ano pedido."""
    hoje = dt.date.today()
    if ano == "atual":
        alvo, janelas = hoje, [(hoje - dt.timedelta(days=200), hoje)]
    else:
        a = int(ano)
        alvo = MARCO if a == 2008 else dt.date(a, 7, 31)
        janelas = [(dt.date(a, 6, 1), dt.date(a, 9, 30)), (dt.date(a, 4, 1), dt.date(a, 11, 30)), (dt.date(a, 1, 1), dt.date(a, 12, 31))]
        if a == 2008:  # marco do Código Florestal: prioriza perto de 22/07/2008
            janelas = [(dt.date(2008, 6, 15), dt.date(2008, 8, 31))] + janelas
    for ini, fim in janelas:
        cands = []
        for prod in PRODUTOS:
            if prod["fim"] and ini > prod["fim"] or fim < prod["ini"]:
                continue
            if ano == "atual" and prod["ms"] != "CB4A-WPM-L4-DN-1":
                continue
            ms = por_data(stac(prod["ms"], bb_folga if not prod["l4"] else bb, ini, fim))
            pan = por_data(stac(prod["pan"][0], bb_folga if not prod["l4"] else bb, ini, fim)) if prod["pan"] else {}
            for d, its in ms.items():
                if prod["pan"] and d not in pan:
                    continue
                if not prod["pan"] and prod["nome"] == "CBERS-2B CCD" and d in por_data(stac("CB2B-HRC-L2-DN-1", bb_folga, d, d)):
                    continue  # esse dia tem HRC: fica com o produto fusionado
                cands.append((prod, d, its, pan.get(d, [])))
        if not cands:
            continue
        # ordem de preferência: resolução, depois proximidade da data alvo
        res = lambda c: c[0]["res_pan"] if c[0]["pan"] else c[0]["res_ms"]
        dist = lambda c: abs((dt.date.fromisoformat(c[1]) - alvo).days)
        cands.sort(key=lambda c: (dist(c) / 45 + res(c) / 4))
        melhor = None
        for c in cands[:10]:
            n = min(nota_nuvem(it, c[0], bb) for it in c[2])
            total = n + dist(c) / 45 + res(c) / 4
            log(f"   {c[1]} {c[0]['nome']}: nuvem/falha {n:.1f}")
            if melhor is None or total < melhor[-1]:
                melhor = (*c, n, total)
            if n < 1.5:
                break
        if melhor and melhor[4] < 25:
            return melhor[:5]
    return None


# ---------------------------------------------------------------- georreferenciamento
def bordas(a):
    a = cv2.GaussianBlur(a, (0, 0), 1.0)
    g = np.hypot(cv2.Sobel(a, cv2.CV_32F, 1, 0), cv2.Sobel(a, cv2.CV_32F, 0, 1))
    m = a > 0
    if m.sum() < 100:
        return g * 0
    return ((g - g[m].mean()) / (g[m].std() + 1e-9) * m).astype(np.float32)


def registrar(ref_pecas, alvo_pecas, crs, centro, meia, res):
    """Desloca 'alvo' para casar com 'ref' (correlação de fase sobre bordas). Devolve (dx, dy, resposta)."""
    x0, y0 = centro[0] - meia, centro[1] + meia
    n = int(2 * meia / res)
    grade = (crs, from_origin(x0, y0, res, res), n, n)
    A = mosaico(ref_pecas, grade, Resampling.average)
    B = mosaico(alvo_pecas, grade, Resampling.average)
    m = (A > 0) & (B > 0)
    if m.mean() < 0.15:
        return 0.0, 0.0, 0.0
    ea, eb = bordas(np.where(m, A, 0)), bordas(np.where(m, B, 0))
    win = cv2.createHanningWindow((n, n), cv2.CV_32F)
    (sx, sy), resp = cv2.phaseCorrelate(ea, eb, win)
    return sx * res, -sy * res, resp


def corrigir(ref, pecas, crs, centro, meia_ampla):
    """Duas passadas: ampla (até ~6 km) a 20 m e fina a 8 m. Devolve (pecas corrigidas, dx, dy, resposta)."""
    dx, dy, r1 = registrar(ref, pecas, crs, centro, meia_ampla, 20)
    p1 = [deslocar(p, dx, dy) for p in pecas]
    dx2, dy2, r2 = registrar(ref, p1, crs, centro, min(meia_ampla, 5000), 8)
    if r2 > 0.05 and math.hypot(dx2, dy2) < 120:
        dx, dy = dx + dx2, dy + dy2
    return [deslocar(p, dx, dy) for p in pecas], dx, dy, max(r1, r2)


# ---------------------------------------------------------------- fusão e realce
def fusao_hpf(pan, ms, razao):
    valido = (pan > 0) & np.all([m > 0 for m in ms], axis=0)
    hp = pan - cv2.GaussianBlur(pan, (0, 0), max(razao / 2.0, 0.6))
    sp = np.std(pan[valido]) if valido.any() else 1
    return [m + hp * (np.std(m[valido]) / (sp + 1e-9)) for m in ms], valido


def realce(bandas, valido, mascara):
    sel = valido & mascara if (valido & mascara).sum() > 500 else valido
    out = []
    for b in bandas:
        lo, hi = np.percentile(b[sel], [1.5, 98.5])
        x = np.clip((b - lo) / (hi - lo + 1e-9), 0, 1) ** 0.9
        out.append(((x * 254 + 1) * valido).astype(np.uint8))
    return out


# ---------------------------------------------------------------- desenho
def fonte(tam, negrito=False):
    for f in (["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"] if negrito else []) + ["/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"]:
        if os.path.exists(f):
            return ImageFont.truetype(f, tam)
    return ImageFont.load_default()


def para_px(geom_utm, tr):
    inv = ~tr
    polis = [geom_utm] if geom_utm.geom_type == "Polygon" else list(getattr(geom_utm, "geoms", []))
    out = []
    for p in polis:
        if p.geom_type != "Polygon":
            continue
        for anel in [p.exterior, *p.interiors]:
            out.append([inv * (x, y) for x, y in anel.coords])
    return out


def tracejado(d, pts, cor, larg, traco=9, vazio=6):
    """Linha tracejada contínua ao longo do contorno (o padrão não reinicia a cada vértice)."""
    ciclo, fase = traco + vazio, 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        L = math.hypot(x1 - x0, y1 - y0)
        p = 0.0
        while p < L:
            pos = (fase + p) % ciclo
            if pos < traco:
                q = min(L, p + traco - pos)
                d.line([(x0 + (x1 - x0) * p / L, y0 + (y1 - y0) * p / L), (x0 + (x1 - x0) * q / L, y0 + (y1 - y0) * q / L)], fill=cor, width=larg)
            else:
                q = min(L, p + ciclo - pos)
            p = q if q > p else p + 0.5
        fase = (fase + L) % ciclo


def escala(d, x, y, mpp, alvo_px, cor=(255, 255, 255)):
    m = alvo_px * mpp
    p = 10 ** math.floor(math.log10(m))
    m = max(v * p for v in (1, 2, 5) if v * p <= m)
    L = m / mpp
    d.rectangle([x, y, x + L, y + 5], fill=cor, outline=(0, 0, 0))
    txt = f"{m/1000:g} km" if m >= 1000 else f"{m:g} m"
    d.text((x + L / 2, y - 3), txt, fill=cor, font=fonte(14, True), anchor="md", stroke_width=2, stroke_fill=(0, 0, 0))


def norte(d, x, y):
    d.polygon([(x, y), (x - 9, y + 24), (x, y + 18)], fill=(255, 255, 255), outline=(0, 0, 0))
    d.polygon([(x, y), (x + 9, y + 24), (x, y + 18)], fill=(40, 40, 40), outline=(0, 0, 0))
    d.text((x, y + 40), "N", fill=(255, 255, 255), font=fonte(15, True), anchor="mm", stroke_width=2, stroke_fill=(0, 0, 0))


# ---------------------------------------------------------------- principal
def main():
    pedido = re.sub(r"[^A-Za-z0-9_-]", "", os.environ.get("PEDIDO", "teste"))[:40] or "teste"
    dados = json.loads(os.environ["DADOS"])
    anos = [a.strip().lower() for a in re.split(r"[,; ]+", os.environ.get("ANOS", "2008, atual")) if a.strip()]
    anos = [a for a in anos if a == "atual" or (a.isdigit() and 2003 <= int(a) <= dt.date.today().year)][:8] or ["2008", "atual"]
    if "atual" not in anos:
        anos.append("atual")  # a cena atual é a referência de georreferenciamento
    imovel = shape(dados["imovel"]).buffer(0)
    alertas = [dict(a, g=shape(a["geom"]).buffer(0)) for a in dados.get("alertas", [])]
    lon, lat = imovel.centroid.x, imovel.centroid.y
    zona = int((lon + 180) // 6) + 1
    epsg = (31960 if lat < 0 else 31954) + zona  # SIRGAS 2000 / UTM
    crs = f"EPSG:{epsg}"
    T = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform
    im_utm = shp_transform(T, imovel)
    for a in alertas:
        a["u"] = shp_transform(T, a["g"])
    bb = list(imovel.bounds)
    bb_folga = [bb[0] - 0.08, bb[1] - 0.08, bb[2] + 0.08, bb[3] + 0.08]
    x0, y0, x1, y1 = im_utm.bounds
    borda = max(300, 0.06 * max(x1 - x0, y1 - y0))
    ext = (x0 - borda, y0 - borda, x1 + borda, y1 + borda)
    centro = ((ext[0] + ext[2]) / 2, (ext[1] + ext[3]) / 2)
    meia_ampla = max(ext[2] - ext[0], ext[3] - ext[1]) / 2 + 6500
    log(f"Pedido {pedido} · CAR {dados.get('car')} · anos {anos} · {crs}")

    os.makedirs("saida", exist_ok=True)
    zpath = f"saida/{pedido}.zip"
    z = zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED)
    ficha, cenas = [], []

    ordem = ["atual"] + [a for a in anos if a != "atual"]
    ref = None
    for ano in ordem:
        log(f"== {ano}")
        esc = escolher(ano, bb, bb_folga)
        if not esc:
            log("   nenhuma cena utilizável")
            ficha.append(f"{ano}: nenhuma cena do INPE sem nuvem sobre o imóvel.")
            continue
        prod, data, its_ms, its_pan, nota = esc
        log(f"   escolhida {data} {prod['nome']} ({len(its_ms)} cena(s) MS, {len(its_pan)} pan)")
        pad = 0.012 if prod["l4"] else 0.07
        ms_p = {k: [p for p in (ler(it["assets"][k]["href"], bb, pad) for it in its_ms if k in it["assets"]) if p] for k in prod["b"]}
        pan_p = [p for p in (ler(it["assets"][prod["pan"][1]]["href"], bb, pad) for it in its_pan) if p] if prod["pan"] else []
        if not all(ms_p.values()):
            ficha.append(f"{ano}: falha ao ler as bandas de {data} ({prod['nome']}).")
            continue
        obs = "nível L4 (ortorretificado pelo INPE)"
        if ano == "atual":  # referência de georreferenciamento: vermelho da cena atual, em área ampla
            ref = [p for p in (ler(it["assets"][prod["b"][1]]["href"], bb, 0.07) for it in its_ms) if p]
        elif not prod["l4"] and not ref:
            obs = "nível L2 SEM correção (não houve cena atual de referência) — conferir o encaixe"
        elif not prod["l4"] and ref:
            ms_corr, dx, dy, r = corrigir(ref, ms_p[prod["b"][1]], crs, centro, meia_ampla)
            ms_p = {k: [deslocar(p, dx, dy) for p in v] for k, v in ms_p.items()}
            obs = f"nível L2 corrigido por registro sobre a cena CBERS-4A atual (deslocamento {math.hypot(dx, dy):.0f} m)"
            if pan_p:
                pan_p, pdx, pdy, pr = corrigir(ref, pan_p, crs, centro, meia_ampla)
                obs += f"; pancromática deslocada {math.hypot(pdx, pdy):.0f} m"
            if r < 0.04:
                obs += " — ATENÇÃO: correlação fraca, conferir o encaixe"
            log("   ", obs)

        res = prod["res_pan"] if pan_p else prod["res_ms"]
        W, H = int(math.ceil((ext[2] - ext[0]) / res)), int(math.ceil((ext[3] - ext[1]) / res))
        grade = (crs, from_origin(ext[0], ext[3], res, res), W, H)
        verde, verm, nir = (mosaico(ms_p[k], grade) for k in prod["b"])
        from rasterio.features import geometry_mask
        mascara = ~geometry_mask([mapping(im_utm)], (H, W), grade[1])
        if pan_p:
            pan = mosaico(pan_p, grade)
            if (pan[mascara] > 0).mean() < 0.9:  # pancromática não cobre o imóvel todo: fica só a multiespectral
                obs += "; pancromática não cobre todo o imóvel — entregue sem fusão"
                pan_p = []
                res = prod["res_ms"]
                W, H = int(math.ceil((ext[2] - ext[0]) / res)), int(math.ceil((ext[3] - ext[1]) / res))
                grade = (crs, from_origin(ext[0], ext[3], res, res), W, H)
                verde, verm, nir = (mosaico(ms_p[k], grade) for k in prod["b"])
                mascara = ~geometry_mask([mapping(im_utm)], (H, W), grade[1])
        if pan_p:
            (verde_f, verm_f, nir_f), valido = fusao_hpf(pan, [verde, verm, nir], prod["res_ms"] / res)
        else:
            verde_f, verm_f, nir_f = verde, verm, nir
            valido = (verde > 0) & (verm > 0) & (nir > 0)
        if valido.mean() < 0.3:
            ficha.append(f"{ano}: a cena {data} cobre só parte do imóvel; descartada.")
            continue
        R, G, B = realce([verm_f, nir_f, verde_f], valido, mascara)  # R3 G4 B2: vermelho, NIR, verde
        rgb = np.stack([R, G, B])
        nome_base = f"{data}_{prod['nome'].split()[0].replace('-', '')}_{('fusao_' if pan_p else '')}{str(res).replace('.', '_')}m"
        # GeoTIFF fusionado
        buf = io.BytesIO()
        with rasterio.MemoryFile() as mf:
            with mf.open(driver="GTiff", width=W, height=H, count=3, dtype="uint8", crs=crs, transform=grade[1],
                         compress="deflate", predictor=2, tiled=True, nodata=0, photometric="RGB") as d:
                d.write(rgb)
                d.update_tags(IMOVEL=dados.get("nome") or "", CAR=dados.get("car") or "", MUNICIPIO=dados.get("mun") or "",
                              SATELITE=prod["nome"], SENSORES=prod["sensores"], DATA=data,
                              CENAS="; ".join(i["id"] for i in its_ms + its_pan), COMPOSICAO="R3 G4 B2 (vermelho, NIR, verde)",
                              FUSAO="HPF (filtro passa-alta), mesma passagem" if pan_p else "sem fusão (sensor sem pancromática)",
                              GEORREF=obs, RESOLUCAO=f"{res} m", SRC=f"SIRGAS 2000 / UTM {zona}{'S' if lat < 0 else 'N'} ({crs})",
                              FONTE="INPE/DGI - Brazil Data Cube", ELABORACAO="R.A AgroTec")
            buf.write(mf.read())
        pasta = "1_imagens_fusionadas" if pan_p else "1_imagens_composicao"
        z.writestr(f"{pasta}/{nome_base}_R3G4B2.tif", buf.getvalue())
        # bandas originais (recortadas, já corrigidas), 16 bits
        for rot, arr in [("B_verde", verde), ("B_vermelho", verm), ("B_NIR", nir)] + ([("PAN_pancromatica", pan)] if pan_p else []):
            b2 = io.BytesIO()
            with rasterio.MemoryFile() as mf:
                with mf.open(driver="GTiff", width=W, height=H, count=1, dtype="uint16", crs=crs, transform=grade[1],
                             compress="deflate", predictor=2, tiled=True, nodata=0) as d:
                    d.write(np.clip(arr, 0, 65535).astype(np.uint16), 1)
                b2.write(mf.read())
            z.writestr(f"2_bandas_originais/{data}_{prod['nome'].split()[0]}/{rot}.tif", b2.getvalue())
        img = Image.fromarray(np.dstack([R, G, B]))
        cenas.append(dict(ano=ano, data=data, prod=prod, res=res, img=img, tr=grade[1], fusao=bool(pan_p)))
        ficha.append(f"{ano}: {dt.date.fromisoformat(data).strftime('%d/%m/%Y')} · {prod['nome']} · {prod['sensores']} · "
                     f"{'fusão HPF a' if pan_p else 'composição a'} {res:g} m · {obs} · cenas: {', '.join(i['id'] for i in its_ms + its_pan)}")
        log("   ok")

    if not cenas:
        z.writestr("LEIA-ME.txt", "Nenhuma cena do INPE utilizável foi encontrada.\n" + "\n".join(ficha))
        z.close()
        return
    cenas.sort(key=lambda c: c["data"])

    # ---------------- mapa da fazenda (cena mais recente)
    titulo_car = f"CAR {dados.get('car') or ''} – {dados.get('mun') or ''}"
    for c in [cenas[-1]] + ([cenas[0]] if len(cenas) > 1 else []):
        im = c["img"].copy()
        esc_ = min(1.0, 2400 / max(im.size))
        if esc_ < 1:
            im = im.resize((int(im.width * esc_), int(im.height * esc_)), Image.LANCZOS)
        tr = c["tr"] * Affine.scale(1 / esc_)
        d = ImageDraw.Draw(im)
        for a in alertas:
            for pts in para_px(a["u"], tr):
                d.line(pts, fill=(255, 40, 40), width=3)
        for pts in para_px(im_utm, tr):
            d.line(pts, fill=(0, 0, 0), width=8)
            d.line(pts, fill=(255, 220, 0), width=5)
        for a in alertas:
            p = a["u"].representative_point()
            x, y = ~tr * (p.x, p.y)
            d.text((x, y), f"A{a['id']}", fill=(255, 255, 255), font=fonte(20, True), anchor="mm", stroke_width=3, stroke_fill=(0, 0, 0))
        escala(d, 30, im.height - 30, c["res"] / esc_, im.width * 0.2)
        norte(d, im.width - 40, 20)
        cab = 70
        fim = Image.new("RGB", (im.width, im.height + cab + 40), "white")
        fim.paste(im, (0, cab))
        d = ImageDraw.Draw(fim)
        d.text((14, 12), f"{titulo_car}  |  Mapa do imóvel", fill=(20, 20, 20), font=fonte(26, True))
        d.text((14, 44), f"{dt.date.fromisoformat(c['data']).strftime('%d/%m/%Y')} – {c['prod']['nome']} {c['res']:g} m  ·  "
               f"Limite do imóvel (amarelo) · Alertas (vermelho)", fill=(60, 60, 60), font=fonte(18))
        d.text((14, fim.height - 30), "Composição R3G4B2 (vegetação em verde, solo/pasto em rosa)  |  Fonte: INPE/DGI · SICAR · PRODES  |  R.A AgroTec",
               fill=(60, 60, 60), font=fonte(16))
        b = io.BytesIO()
        fim.save(b, "PNG", optimize=True)
        z.writestr(f"3_mapa_da_fazenda_{c['data']}.png", b.getvalue())

    # ---------------- pranchas de zoom por alerta (padrão R.A AgroTec)
    P, cols, gap, cab = 300, 4, 9, 84
    for a in alertas:
        ax0, ay0, ax1, ay1 = a["u"].bounds
        lado = max(ax1 - ax0, ay1 - ay0) * 1.7 + 200
        cx, cy = (ax0 + ax1) / 2, (ay0 + ay1) / 2
        n = len(cenas)
        rows = math.ceil(n / cols)
        Wp = cols * P + (cols - 1) * gap + 14
        Hp = cab + rows * (P + 34) + 40
        folha = Image.new("RGB", (Wp, Hp), "white")
        d = ImageDraw.Draw(folha)
        d.text((7, 10), f"Alerta A{a['id']} – {a.get('fonte') or 'PRODES'}{(' d' + str(a['ano'])) if a.get('ano') else ''} – {str(round(a.get('ha') or 0, 1)).replace('.', ',')} ha",
               fill=(20, 20, 20), font=fonte(20, True))
        d.text((7, 38), titulo_car, fill=(70, 70, 70), font=fonte(15))
        for i, c in enumerate(cenas):
            col, row = i % cols, i // cols
            X, Y = 7 + col * (P + gap), cab + row * (P + 34)
            res = c["res"]
            inv = ~c["tr"]
            px0, py0 = inv * (cx - lado / 2, cy + lado / 2)
            px1, py1 = inv * (cx + lado / 2, cy - lado / 2)
            rec = c["img"].crop((int(px0), int(py0), int(px1), int(py1))).resize((P, P), Image.BICUBIC)
            k = P / (px1 - px0)
            dd = ImageDraw.Draw(rec)
            trz = Affine.translation(cx - lado / 2, cy + lado / 2) * Affine.scale(lado / P, -lado / P)
            for pts in para_px(a["u"], trz):
                tracejado(dd, pts, (0, 0, 0), 5)
                tracejado(dd, pts, (255, 255, 255), 3)
            escala(dd, 10, 22, lado / P, P * 0.15)
            folha.paste(rec, (X, Y + 22))
            d.text((X, Y), f"{dt.date.fromisoformat(c['data']).strftime('%d/%m/%Y')} – {c['prod']['nome'].split()[0]} {res:g} m",
                   fill=(20, 20, 20), font=fonte(14, True))
        y = Hp - 26
        tracejado(d, [(10, y + 8), (34, y + 8)], (0, 0, 0), 3, 5, 3)
        d.text((42, y), "Polígono do alerta", fill=(40, 40, 40), font=fonte(13))
        d.text((Wp - 8, y), "Composição R3G4B2 (vegetação em verde, solo/pasto em rosa)  |  Fonte: INPE/DGI", fill=(40, 40, 40), font=fonte(13), anchor="ra")
        b = io.BytesIO()
        folha.save(b, "PNG", optimize=True)
        z.writestr(f"4_pranchas_por_alerta/Alerta_A{a['id']}.png", b.getvalue())

    # ---------------- vetores e ficha técnica
    z.writestr("5_limite_e_alertas.geojson", json.dumps({"type": "FeatureCollection", "features":
        [{"type": "Feature", "properties": {"tipo": "imovel", "car": dados.get("car")}, "geometry": mapping(imovel)}] +
        [{"type": "Feature", "properties": {"tipo": "alerta", "id": f"A{a['id']}", "fonte": a.get("fonte"), "ano": a.get("ano"), "ha": a.get("ha")},
          "geometry": mapping(a["g"])} for a in alertas]}, ensure_ascii=False))
    z.writestr("FICHA_TECNICA.txt", "\r\n".join([
        "R.A AgroTec – Soluções Rurais · Imagens do INPE georreferenciadas",
        f"Imóvel: {dados.get('nome') or '—'}", f"CAR: {dados.get('car') or '—'}", f"Município: {dados.get('mun') or '—'}",
        f"Pedido: {pedido} · gerado em {dt.date.today().strftime('%d/%m/%Y')}", "",
        f"Sistema de referência: SIRGAS 2000 / UTM zona {zona}{'S' if lat < 0 else 'N'} ({crs})",
        "Composição: R3 G4 B2 — R = vermelho (0,63–0,69 µm), G = infravermelho próximo (0,77–0,89 µm), B = verde (0,52–0,59 µm).",
        "Vegetação aparece em verde; solo exposto, pasto e palha em rosa/magenta; água em azul escuro.",
        "Fusão: HPF (filtro passa-alta) — detalhe da banda pancromática injetado nas bandas multiespectrais da mesma passagem.",
        "Georreferenciamento: produtos L4 do INPE são ortorretificados; produtos L2 (CBERS-2/2B e LISS-III) foram corrigidos",
        "por correlação de fase sobre a cena CBERS-4A L4 mais recente.",
        "Fonte das imagens: INPE/DGI – Brazil Data Cube (data.inpe.br). Limite: SICAR. Alertas: PRODES (INPE) e MapBiomas Alerta.", "",
        "CENAS:", *ficha, "",
        "PASTAS:",
        "1_imagens_fusionadas / 1_imagens_composicao: GeoTIFF RGB 8 bits, pronto para QGIS/ArcGIS/Google Earth Pro.",
        "2_bandas_originais: bandas recortadas em 16 bits (verde, vermelho, NIR e pancromática), já corrigidas.",
        "3_mapa_da_fazenda: imagem com limite do imóvel e alertas numerados.",
        "4_pranchas_por_alerta: cenas históricas em zoom sobre cada alerta.",
        "5_limite_e_alertas.geojson: limite do imóvel e alertas.",
    ]))
    z.close()
    log("Pronto:", zpath, os.path.getsize(zpath) // 1024, "KB")


if __name__ == "__main__":
    main()
