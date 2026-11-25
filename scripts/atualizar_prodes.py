#!/usr/bin/env python3
"""Baixa os polígonos do PRODES (INPE/TerraBrasilis) e guarda no site em blocos de 0,25°.

O app lê só os blocos que cobrem o imóvel, sem depender do servidor do INPE.
Uma camada só é regravada quando o conteúdo mudou; se o download falhar,
a cópia anterior continua valendo.

Uso: python3 scripts/atualizar_prodes.py [camada ...]
"""
import hashlib
import json
import math
import os
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date

TB = "https://terrabrasilis.dpi.inpe.br/geoserver/ows"
RAIZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "prodes")
PASSO = 0.25          # tamanho do bloco em graus
PAGINA = 5000         # polígonos por requisição
ANO_MIN = 2008        # antes de 22/07/2008 é área rural consolidada: não precisa guardar
UA = "Mozilla/5.0 (compatible; RA-AgroTec-atualizacao/1.0)"

CAMADAS = {
    "cerrado": "prodes-cerrado-nb:yearly_deforestation",
    "legal_amz": "prodes-legal-amz:yearly_deforestation",
    "amazonia_bioma": "prodes-amazon-nb:yearly_deforestation_biome",
    "mata_atlantica": "prodes-mata-atlantica-nb:yearly_deforestation",
    "caatinga": "prodes-caatinga-nb:yearly_deforestation",
    "pantanal": "prodes-pantanal-nb:yearly_deforestation",
    "pampa": "prodes-pampa-nb:yearly_deforestation",
}
# campos usados pelo app (o resto é descartado para economizar espaço)
CAMPOS = ["year", "ano", "image_date", "date", "view_date", "uid", "fid", "gid", "id",
          "main_class", "class_name", "state", "area_km", "satellite", "sensor", "path_row"]


def baixar(url, tentativas=8):
    espera = 20
    for t in range(tentativas):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.loads(r.read())
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ConnectionError) as e:
            print(f"    falha ({e}); nova tentativa em {espera}s", flush=True)
            time.sleep(espera)
            espera = min(espera * 2, 600)
    raise RuntimeError("servidor do INPE não respondeu")


def arred(c):
    if isinstance(c[0], (int, float)):
        return [round(c[0], 5), round(c[1], 5)]
    out = [arred(x) for x in c]
    if out and isinstance(out[0][0], float):  # anel: tira pontos repetidos
        anel = [out[0]]
        for p in out[1:]:
            if p != anel[-1]:
                anel.append(p)
        return anel
    return out


def bbox(c, b=None):
    b = b or [999, 999, -999, -999]
    if isinstance(c[0], (int, float)):
        b[0], b[1], b[2], b[3] = min(b[0], c[0]), min(b[1], c[1]), max(b[2], c[0]), max(b[3], c[1])
    else:
        for x in c:
            bbox(x, b)
    return b


def atualizar(nome, camada):
    print(f"== {nome} ({camada})", flush=True)
    blocos, n, anos, ini = {}, 0, set(), 0
    while True:
        q = urllib.parse.urlencode({
            "service": "WFS", "version": "2.0.0", "request": "GetFeature", "typeNames": camada,
            "outputFormat": "application/json", "srsName": "EPSG:4326",
            "count": PAGINA, "startIndex": ini,
        })
        d = baixar(f"{TB}?{q}")
        fs = d.get("features") or []
        for f in fs:
            g = f.get("geometry")
            if not g or not g.get("coordinates"):
                continue
            p = {k: v for k, v in (f.get("properties") or {}).items() if k in CAMPOS and v not in (None, "")}
            p["_id"] = f"{nome}.{f.get('id') or n}"
            a = p.get("year") or p.get("ano")
            if a:
                try:
                    a = int(float(a))
                except ValueError:
                    a = None
            if a and a < ANO_MIN:
                continue
            if a:
                anos.add(a)
            g = {"type": g["type"], "coordinates": arred(g["coordinates"])}
            b = bbox(g["coordinates"])
            feat = {"type": "Feature", "properties": p, "geometry": g}
            for tx in range(math.floor(b[0] / PASSO), math.floor(b[2] / PASSO) + 1):
                for ty in range(math.floor(b[1] / PASSO), math.floor(b[3] / PASSO) + 1):
                    blocos.setdefault(f"{tx}_{ty}", []).append(feat)
            n += 1
        print(f"    página {ini // PAGINA + 1}: {n} polígonos de {ANO_MIN} em diante", flush=True)
        if len(fs) < PAGINA:
            break
        ini += PAGINA
        time.sleep(2)  # não sobrecarregar o INPE
    if n == 0:
        raise RuntimeError("camada veio vazia")

    # monta em pasta temporária e só troca se o conteúdo mudou
    h = hashlib.sha256()
    for k in sorted(blocos):
        h.update(k.encode())
        h.update(json.dumps(blocos[k], separators=(",", ":"), sort_keys=True).encode())
    assinatura = h.hexdigest()[:16]
    idx = ler_indice()
    if idx.get(nome, {}).get("assinatura") == assinatura:
        idx[nome]["conferido"] = date.today().isoformat()
        gravar_indice(idx)
        print("    sem mudanças", flush=True)
        return False
    tmp = os.path.join(RAIZ, f".{nome}.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    for k, fs in blocos.items():
        with open(os.path.join(tmp, f"{k}.json"), "w") as fh:
            json.dump({"type": "FeatureCollection", "features": fs}, fh, separators=(",", ":"))
    dest = os.path.join(RAIZ, nome)
    shutil.rmtree(dest, ignore_errors=True)
    os.rename(tmp, dest)
    idx[nome] = {"camada": camada, "versao": date.today().isoformat(), "poligonos": n,
                 "anos": [min(anos), max(anos)] if anos else None, "blocos": sorted(blocos), "ano_min": ANO_MIN,
                 "passo": PASSO, "assinatura": assinatura, "conferido": date.today().isoformat()}
    gravar_indice(idx)
    print(f"    gravado: {n} polígonos em {len(blocos)} blocos", flush=True)
    return True


def ler_indice():
    try:
        with open(os.path.join(RAIZ, "index.json")) as fh:
            return json.load(fh)
    except FileNotFoundError:
        return {}


def gravar_indice(idx):
    os.makedirs(RAIZ, exist_ok=True)
    with open(os.path.join(RAIZ, "index.json"), "w") as fh:
        json.dump(idx, fh, separators=(",", ":"), ensure_ascii=False)


def main():
    nomes = sys.argv[1:] or list(CAMADAS)
    os.makedirs(RAIZ, exist_ok=True)
    # roda todo dia, mas só baixa tudo de novo se faltar camada ou se a base tiver mais de 25 dias
    if not os.environ.get("FORCAR"):
        idx = ler_indice()
        faltam = [n for n in nomes if n not in idx]
        if not faltam:
            idade = min((date.today() - date.fromisoformat(idx[n]["conferido"] if "conferido" in idx[n] else idx[n]["versao"])).days for n in nomes)
            if idade < 25:
                print(f"Base conferida há {idade} dia(s); nada a fazer.")
                return
        else:
            nomes = faltam + [n for n in nomes if n not in faltam]
    falhas = []
    for nome in nomes:
        try:
            atualizar(nome, CAMADAS[nome])
        except Exception as e:  # uma camada com problema não derruba as outras
            print(f"    ERRO em {nome}: {e} — mantida a cópia anterior", flush=True)
            falhas.append(nome)
    if falhas:
        print("Camadas não atualizadas:", ", ".join(falhas))
    if len(falhas) == len(nomes):
        sys.exit(1)


if __name__ == "__main__":
    main()
