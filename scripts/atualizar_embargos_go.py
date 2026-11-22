"""Atualiza a cópia dos embargos estaduais de Goiás (SEMAD/SIGA-GO). Sem nome/CPF."""
import json,urllib.request,time,sys,os
UA='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36'
B="https://siga.meioambiente.go.gov.br/geoserver/ows?service=WFS&version=2.0.0&request=GetFeature&typeNames=geonode:vw_embargo_ina_sga_poligono_publico&srsName=EPSG:4326&outputFormat=application/json&sortBy=codigo&count=3000&startIndex="
feats=[];i=0
while True:
    for t in range(4):
        try:
            d=json.load(urllib.request.urlopen(urllib.request.Request(B+str(i),headers={'User-Agent':UA}),timeout=180));break
        except Exception as e:
            print('retry',i,e,file=sys.stderr);time.sleep(5)
    fs=d.get('features',[]);feats+=fs;print(i,len(fs),len(feats),file=sys.stderr)
    if len(fs)<3000:break
    i+=3000
json.dump({'type':'FeatureCollection','features':feats},open('/tmp/bruto.json','w'))
print(len(feats))

# ---- processamento ----
import json,math,os,re,collections,sys
from shapely.geometry import shape,mapping
from shapely.validation import make_valid
from shapely.ops import unary_union
from datetime import date
GO=(-53.6,-19.8,-45.4,-12.2)  # extensão de Goiás com margem
def dentro(b):return b[0]>=GO[0] and b[2]<=GO[2] and b[1]>=GO[1] and b[3]<=GO[3]
def clean(v):return None if v in (None,'None','') else str(v).strip()
def rnd(o):return [rnd(x) for x in o] if isinstance(o,(list,tuple)) else round(o,6)
d=json.load(open('/tmp/bruto.json'))
tiles=collections.defaultdict(list);n=0;fora=0;troca=0;ruim=0;S=0.5
for f in d['features']:
    if not f.get('geometry'):ruim+=1;continue
    try:g=shape(f['geometry'])
    except Exception:ruim+=1;continue
    if g.is_empty:ruim+=1;continue
    b=g.bounds
    if not dentro(b):
        from shapely.ops import transform
        g2=transform(lambda x,y,z=None:(y,x),g)
        if dentro(g2.bounds):g=g2;troca+=1
        else:fora+=1;continue
    try:
        if not g.is_valid:g=make_valid(g)
        if g.geom_type not in('Polygon','MultiPolygon'):
            ps=[p for p in getattr(g,'geoms',[]) if p.geom_type in('Polygon','MultiPolygon')]
            if not ps:ruim+=1;continue
            g=unary_union(ps)
        g=g.simplify(0.00002,preserve_topology=True)
        if g.is_empty:ruim+=1;continue
    except Exception:ruim+=1;continue
    p=f['properties'];ns=p.get('numero_processo_sei_sga') or ''
    sei=re.search(r'SEI = (\d+)',ns);sga=re.search(r'SGA = (\d+)',ns)
    props={k:v for k,v in {'codigo':clean(p.get('codigo')),'data':(clean(p.get('data_criacao')) or '')[:10] or None,'sei':sei.group(1) if sei else None,
      'sga':sga.group(1) if sga else clean(p.get('codigo_processo')),'categoria':clean(p.get('categoria')),'sub':clean(p.get('sub_categoria')),
      'municipio':clean(p.get('municipio_nome')),'area':clean(p.get('area')),'origem':clean(p.get('origem'))}.items() if v}
    gm=mapping(g);ft={'type':'Feature','properties':props,'geometry':{'type':gm['type'],'coordinates':rnd(gm['coordinates'])}};n+=1
    x0,y0,x1,y1=g.bounds
    for tx in range(math.floor(x0/S),math.floor(x1/S)+1):
        for ty in range(math.floor(y0/S),math.floor(y1/S)+1):tiles[(tx,ty)].append(ft)
import shutil
shutil.rmtree('emb_go',ignore_errors=True);os.makedirs('emb_go',exist_ok=True)
for (tx,ty),fs in tiles.items():
    json.dump({'type':'FeatureCollection','features':fs},open(f'emb_go/{tx}_{ty}.json','w'),separators=(',',':'),ensure_ascii=False)
json.dump({'fonte':'SEMAD-GO / SIGA-GO — Embargos Estaduais polígono (INÃ/SGA, público)','camada':'geonode:vw_embargo_ina_sga_poligono_publico','atualizado':date.today().isoformat(),'total':n,'grade':S,'tiles':sorted(f'{a}_{b}' for a,b in tiles)},open('emb_go/index.json','w'),ensure_ascii=False)
sz=[os.path.getsize('emb_go/'+x) for x in os.listdir('emb_go')]
print(f'feições {n} | eixos trocados corrigidos {troca} | fora de GO descartadas {fora} | inválidas {ruim} | blocos {len(tiles)} | total {sum(sz)/1e6:.1f} MB | maior {max(sz)/1e6:.2f} MB')
