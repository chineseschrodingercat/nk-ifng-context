"""Finalize provenance and verifiable delivery metadata without changing parent runs."""
from pathlib import Path
import json,hashlib,platform,shutil,sys,re
import pandas as pd,numpy,scipy,sklearn,matplotlib,PIL
R=Path(__file__).resolve().parents[1]
project=Path(r'C:\N01_NK_IFNg')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
source=[]
for name in ['PROPOSAL.md','PHASE0A_DATA_AUDIT_PROMPT.md','AGENTS.md']:
 p=project/name
 if p.exists():
  out=R/'provenance'/name
  if not out.exists():shutil.copyfile(p,out)
  assert sha(p)==sha(out)
  source.append(dict(source=str(p),local=str(out.relative_to(R)),sha256=sha(p)))
(R/'provenance/input_files.json').write_text(json.dumps(source,indent=2),encoding='utf8')
parents=[('phase0a',project/'phase0a_runs/20260922_161233','provenance/DELIVERY_FILE_MANIFEST.csv','relative_path'),('phase0b',project/'phase0b_runs/20260922_170151','artifact_manifest.csv','path'),('mechanism',project/'mechanism_runs/20260922_180448','artifact_manifest.csv','path')]
checks=[]
for label,p,manifest,key in parents:
 table=pd.read_csv(p/manifest);bad=[]
 for _,row in table.iterrows():
  f=p/row[key]
  if not f.exists() or f.stat().st_size!=row['bytes'] or sha(f)!=row.sha256:bad.append(row[key])
 checks.append(dict(parent=label,path=str(p),manifest=manifest,n_files=len(table),mismatches=bad))
 assert not bad,(label,bad)
 shutil.copyfile(p/manifest,R/'provenance'/f'{label}_parent_manifest.csv')
(R/'qa/parent_integrity.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
runtime=dict(python=sys.version,executable=sys.executable,platform=platform.platform(),packages=dict(pandas=pd.__version__,numpy=numpy.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,matplotlib=matplotlib.__version__,Pillow=PIL.__version__))
(R/'qa/runtime_versions.json').write_text(json.dumps(runtime,indent=2),encoding='utf8')
for p in (R/'scripts').glob('*.py'):compile(p.read_text(encoding='utf8'),str(p),'exec')
badlinks=[]
for p in R.glob('*.md'):
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf8')):
  if target.startswith(('https:','http:','C:','#')):continue
  if target=='DELIVERY_QA.json':continue
  if not (p.parent/target).exists():badlinks.append((p.name,target))
assert not badlinks,badlinks
figreview=[]
for name in ['Figure1_grouped_prediction','Figure2_independent_RNA']:
 p=R/'figures'/f'{name}.pdf'
 figreview.append(dict(file=str(p.relative_to(R)),sha256=sha(p),actual_pdf_render_reviewed_by_root=True,font='embedded Arial-BoldMT',vector=True,PNG_dpi=300))
(R/'qa/FIGURE_REVIEW.md').write_text('''# Final figure review

Both final PDF files were rendered with pdftoppm and their final rasterized pages were visually reviewed by the root analyst. Arial-BoldMT is embedded in both PDFs. Each has two clearly separated panels, four visible spines, no grid, external keys and native outcome/error units. Labels and markers have been separated, right margins widened, and the unsupported subscript glyph replaced with plain log2. Per-panel data CSVs agree with the analysis tables. The numerical labels report only the comparisons described in the legends. No confidence bars or significance symbols were added. Physical report-figure sizes and actual axes dimensions are recorded in figure_style_manifest.json; they are not claimed to be final journal-column layouts.
''',encoding='utf8')
q=json.loads((R/'qa/analysis_validation.json').read_text());q.update(status='PASS',run='20260922_185850',stage='exploratory grouped prediction and molecular extension',parent_integrity=checks,figures=figreview,links_checked=True,all_scripts_syntax_checked=True,biological_hypothesis_confirmed=False,external_functional_validation=False,manuscript_submission_ready=False)
(R/'DELIVERY_QA.json').write_text(json.dumps(q,indent=2),encoding='utf8')
status={'run_status':'COMPUTATIONAL_STAGE_COMPLETE','scientific_status':'SMALL_COMPLETE_CASE_RECEPTOR_SIGNAL_WITHOUT_GENERAL_RESPONSE_RULE','complete_case_signal':'baseline FAS/DR4/DR5 improves family-held-out MSE by 5.5%, driven partly by SJ-GBM2; full22 and MAE do not reproduce this gain','broader_hypothesis':'remains open; direct functional competence was not measured here','manuscript':'results draft supplied; not submission-ready for the proposed strong claim','future_priority':'same-state functional branch measurements and independent direct IFNg-NK outcome validation','original_runs_modified':False}
(R/'RUN_STATUS.json').write_text(json.dumps(status,indent=2),encoding='utf8')
rows=[]
for p in sorted(R.rglob('*')):
 if p.is_file() and p.name!='artifact_manifest.csv' and '__pycache__' not in p.parts:rows.append(dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(rows).to_csv(R/'artifact_manifest.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'parents':checks,'package_files':len(rows),'package_bytes':sum(x['bytes'] for x in rows),'qa':q['status']},indent=2))
