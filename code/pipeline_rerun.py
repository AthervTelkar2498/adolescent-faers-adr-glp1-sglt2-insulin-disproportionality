#!/usr/bin/env python3
"""
FAERS Adolescent Pharmacovigilance Pipeline — Memory-Efficient Re-Run
======================================================================
Tasks A–D per ANTIGRAVITY_TASK_BRIEF.md.
Key fix: filter DRUG table to panel ingredients BEFORE upper-casing all rows.
"""

import pandas as pd
import numpy as np
import os, json, warnings
from math import sqrt, exp, log2
from pathlib import Path

warnings.filterwarnings('ignore')

# ===================================================================
# CONFIG
# ===================================================================
DATA_ROOT  = r"e:\Adverse drugs Project 2\faers-pediatric-adr\data"
OUTPUT_DIR = r"e:\Adverse drugs Project 2\faers-pediatric-adr\pipeline_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Drug panels (all uppercase, as they appear in prod_ai / drugname) ---
OBESITY_BROAD    = ['METFORMIN','ATORVASTATIN','LISINOPRIL','LOSARTAN',
                    'SEMAGLUTIDE','LIRAGLUTIDE','TIRZEPATIDE',
                    'DAPAGLIFLOZIN','EMPAGLIFLOZIN','ORLISTAT',
                    'INSULIN GLARGINE','INSULIN LISPRO',
                    'INSULIN ASPART','INSULIN DEGLUDEC']
DIABETES_BROAD   = ['METFORMIN','INSULIN GLARGINE','INSULIN LISPRO',
                    'INSULIN ASPART','INSULIN DEGLUDEC','SEMAGLUTIDE',
                    'DAPAGLIFLOZIN','EMPAGLIFLOZIN','SITAGLIPTIN','PIOGLITAZONE']
OBESITY_SEL      = ['SEMAGLUTIDE','LIRAGLUTIDE','TIRZEPATIDE','ORLISTAT']
DIABETES_SEL     = ['SEMAGLUTIDE','DAPAGLIFLOZIN','EMPAGLIFLOZIN','SITAGLIPTIN']
ALL_PANEL_DRUGS  = set(OBESITY_BROAD + DIABETES_BROAD)   # union

# Partial name fragments to pre-filter DRUG rows cheaply before exact matching
# This avoids upper-casing 39M rows; we only upper-case the small subset
PANEL_FRAGMENTS = [
    'metformin','atorvastatin','lisinopril','losartan',
    'semaglutide','liraglutide','tirzepatide','dapagliflozin',
    'empagliflozin','orlistat','insulin','sitagliptin','pioglitazone'
]

def _fragment_filter(series):
    """Return boolean mask: row contains ANY panel fragment (case-insensitive)."""
    s = series.fillna('').str.lower()
    mask = pd.Series(False, index=series.index)
    for frag in PANEL_FRAGMENTS:
        mask |= s.str.contains(frag, regex=False)
    return mask


# ===================================================================
# FILE DISCOVERY
# ===================================================================
def discover_files(prefix):
    files = []
    for year in ['2021','2022','2023','2024']:
        ydir = os.path.join(DATA_ROOT, year)
        if not os.path.isdir(ydir):
            continue
        for qdir in sorted(os.listdir(ydir)):
            adir = os.path.join(ydir, qdir, 'ASCII')
            if not os.path.isdir(adir):
                continue
            for f in os.listdir(adir):
                if f.upper().startswith(prefix.upper()) and f.endswith('.txt'):
                    files.append(os.path.join(adir, f))
    # 2025 consolidated
    cdir = os.path.join(DATA_ROOT, f'{prefix.upper()}24-25')
    if os.path.isdir(cdir):
        for f in os.listdir(cdir):
            if f.upper().startswith(prefix.upper()) and f.endswith('.txt') and '25Q' in f.upper():
                files.append(os.path.join(cdir, f))
    return sorted(files)


def load_table(files, usecols=None, dtype=str):
    dfs = []
    for fp in files:
        name = os.path.basename(fp)
        try:
            df = pd.read_csv(fp, sep='$', dtype=dtype, usecols=usecols,
                             encoding='latin-1', on_bad_lines='skip', low_memory=False)
            df.columns = [c.lower().strip() for c in df.columns]
            print(f"  {name}: {len(df):,} rows")
            dfs.append(df)
        except Exception as e:
            print(f"  {name}: ERROR {e}")
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


# ===================================================================
# AGE CONVERSION
# ===================================================================
def to_years(age, cod):
    try:
        a = float(age)
        if a <= 0:
            return None
    except:
        return None
    c = str(cod).strip().upper() if pd.notna(cod) else ''
    return {'YR': a, 'MON': a/12, 'WK': a/52, 'DY': a/365.25,
            'HR': a/8766, 'DEC': a*10}.get(c, None)


# ===================================================================
# TASK A
# ===================================================================
def task_a():
    print("="*68)
    print("TASK A: REBUILD COHORT")
    print("="*68)

    # --- DEMO ---
    print("\n[A1] Loading DEMO files...")
    demo_files = discover_files('DEMO')
    print(f"  {len(demo_files)} files found")
    demo = load_table(demo_files,
                      usecols=['primaryid','caseid','caseversion',
                               'age','age_cod','sex','wt','wt_cod',
                               'reporter_country','occr_country','init_fda_dt'])
    print(f"  Total raw rows: {len(demo):,}")

    demo['caseid']      = demo['caseid'].astype(str).str.strip()
    demo['primaryid']   = demo['primaryid'].astype(str).str.strip()
    demo['caseversion'] = pd.to_numeric(demo['caseversion'], errors='coerce').fillna(0).astype(int)

    # De-duplicate: keep latest caseversion per caseid
    demo = demo.sort_values(['caseid','caseversion']).drop_duplicates('caseid', keep='last')
    n_dedup = len(demo)
    print(f"  After dedup: {n_dedup:,} unique cases")

    # Age conversion — no imputation
    demo['age_years'] = demo.apply(lambda r: to_years(r['age'], r['age_cod']), axis=1)
    n_missing = demo['age_years'].isna().sum()
    print(f"  Missing/unknown age dropped: {n_missing:,}")

    adol = demo.dropna(subset=['age_years'])
    adol = adol[(adol['age_years'] >= 12) & (adol['age_years'] < 18)].copy()
    n_adol = len(adol)
    print(f"  Adolescents (12–<18): {n_adol:,}")
    print(f"  Age  mean={adol['age_years'].mean():.2f}  SD={adol['age_years'].std():.2f}")
    print(f"  Female {(adol['sex']=='F').mean()*100:.1f}%  Male {(adol['sex']=='M').mean()*100:.1f}%")

    adol_ids = set(adol['caseid'].values)

    # Save
    adol.to_csv(os.path.join(OUTPUT_DIR,'adolescent_cohort.csv'), index=False)
    stats = {
        'total_raw_demo_rows': int(len(demo) + n_missing),
        'deduplicated_cases': int(n_dedup),
        'dropped_missing_age': int(n_missing),
        'adolescent_cases': int(n_adol),
        'age_mean': round(float(adol['age_years'].mean()),2),
        'age_sd':   round(float(adol['age_years'].std()),2),
        'female_pct': round(float((adol['sex']=='F').mean()*100),1),
    }

    # --- DRUG (memory-efficient) ---
    print("\n[A2] Loading DRUG files (panel-filtered)...")
    drug_files = discover_files('DRUG')
    print(f"  {len(drug_files)} files found")
    drug_dfs = []
    for fp in drug_files:
        name = os.path.basename(fp)
        try:
            df = pd.read_csv(fp, sep='$', dtype=str,
                             usecols=['primaryid','caseid','drug_seq','role_cod','drugname','prod_ai'],
                             encoding='latin-1', on_bad_lines='skip', low_memory=False)
            df.columns = [c.lower().strip() for c in df.columns]
            # Step 1: filter to adolescent cases first (cheap string compare)
            df['caseid'] = df['caseid'].astype(str).str.strip()
            df = df[df['caseid'].isin(adol_ids)]
            if len(df) == 0:
                print(f"  {name}: 0 adolescent rows")
                continue
            # Step 2: filter to panel fragments before upper-casing
            ai_mask = _fragment_filter(df['prod_ai'].fillna('') + ' ' + df['drugname'].fillna(''))
            df = df[ai_mask].copy()
            if len(df) == 0:
                print(f"  {name}: 0 panel rows after fragment filter")
                continue
            # Step 3: now it's small enough to upper-case safely
            df['ingredient'] = (df['prod_ai'].fillna(df['drugname'])).str.upper().str.strip()
            df['ingredient'] = df['ingredient'].str.replace(r'\s+',' ',regex=True)
            df['role_cod']   = df['role_cod'].str.strip().str.upper()
            # Step 4: exact panel match
            df = df[df['ingredient'].isin(ALL_PANEL_DRUGS)].copy()
            print(f"  {name}: {len(df):,} panel rows")
            drug_dfs.append(df)
        except Exception as e:
            print(f"  {name}: ERROR {e}")

    drug_adol = pd.concat(drug_dfs, ignore_index=True) if drug_dfs else pd.DataFrame()
    print(f"  Total panel drug rows: {len(drug_adol):,}")
    drug_adol.to_csv(os.path.join(OUTPUT_DIR,'drug_panel_adolescent.csv'), index=False)

    # --- REAC ---
    print("\n[A3] Loading REAC files...")
    reac_files = discover_files('REAC')
    print(f"  {len(reac_files)} files found")
    reac_dfs = []
    for fp in reac_files:
        name = os.path.basename(fp)
        try:
            df = pd.read_csv(fp, sep='$', dtype=str,
                             usecols=['primaryid','caseid','pt'],
                             encoding='latin-1', on_bad_lines='skip', low_memory=False)
            df.columns = [c.lower().strip() for c in df.columns]
            df['caseid'] = df['caseid'].astype(str).str.strip()
            df = df[df['caseid'].isin(adol_ids)]
            df['pt'] = df['pt'].str.strip()
            reac_dfs.append(df)
            print(f"  {name}: {len(df):,} rows")
        except Exception as e:
            print(f"  {name}: ERROR {e}")

    reac_adol = pd.concat(reac_dfs, ignore_index=True) if reac_dfs else pd.DataFrame()
    print(f"  Total REAC rows: {len(reac_adol):,}")
    reac_adol.to_csv(os.path.join(OUTPUT_DIR,'reac_adolescent.csv'), index=False)

    # --- OUTC ---
    print("\n[A4] Loading OUTC files...")
    outc_files = discover_files('OUTC')
    outc_dfs = []
    for fp in outc_files:
        try:
            df = pd.read_csv(fp, sep='$', dtype=str,
                             usecols=['primaryid','caseid','outc_cod'],
                             encoding='latin-1', on_bad_lines='skip', low_memory=False)
            df.columns = [c.lower().strip() for c in df.columns]
            df['caseid'] = df['caseid'].astype(str).str.strip()
            df = df[df['caseid'].isin(adol_ids)]
            outc_dfs.append(df)
        except Exception as e:
            print(f"  ERROR {e}")
    outc_adol = pd.concat(outc_dfs, ignore_index=True) if outc_dfs else pd.DataFrame()
    print(f"  Total OUTC rows: {len(outc_adol):,}")
    outc_adol.to_csv(os.path.join(OUTPUT_DIR,'outc_adolescent.csv'), index=False)

    # --- (caseid, ingredient, PT) triples — deduplicated ---
    print("\n[A5] Building (caseid, ingredient, PT) triples...")
    if len(drug_adol) > 0 and len(reac_adol) > 0:
        pairs = drug_adol[['caseid','ingredient','role_cod']].merge(
            reac_adol[['caseid','pt']], on='caseid', how='inner')
        pairs = pairs.drop_duplicates(subset=['caseid','ingredient','pt'])
        print(f"  Unique triples: {len(pairs):,}")
        pairs.to_csv(os.path.join(OUTPUT_DIR,'drug_event_pairs.csv'), index=False)
    else:
        pairs = pd.DataFrame()
        print("  No pairs (empty drug or reac table)")

    with open(os.path.join(OUTPUT_DIR,'cohort_stats.json'),'w') as f:
        json.dump(stats, f, indent=2)

    return adol, drug_adol, reac_adol, outc_adol, pairs, stats


# ===================================================================
# TASK B: DISPROPORTIONALITY
# ===================================================================
def ror(a,b,c,d):
    if min(a,b,c,d)==0: return None,None,None
    r = (a*d)/(b*c); se=sqrt(1/a+1/b+1/c+1/d)
    return r, exp(np.log(r)-1.96*se), exp(np.log(r)+1.96*se)

def prr(a,b,c,d):
    if a==0 or (a+b)==0 or c==0 or (c+d)==0: return None,None,None
    p = (a/(a+b))/(c/(c+d)); se=sqrt(1/a-1/(a+b)+1/c-1/(c+d))
    return p, exp(np.log(p)-1.96*se), exp(np.log(p)+1.96*se)

def ic_bcpnn(a,b,c,d):
    N=a+b+c+d
    if N==0 or a==0: return None,None
    E=((a+b)*(a+c))/N
    if E==0: return None,None
    ic=log2(a/E)
    ic025=ic-1.96/(sqrt(a)*np.log(2))
    return ic,ic025

def ebgm(a,b,c,d):
    N=a+b+c+d
    if N==0 or a==0: return None,None
    E=((a+b)*(a+c))/N
    if E==0: return None,None
    eb=(a+0.5)/(E+0.5)
    eb05=eb*exp(-1.645/sqrt(a+0.5))
    return eb,eb05

def run_disprop(adol, drug_adol, reac_adol, pairs, panel_drugs, panel_name, role):
    total = len(adol)
    print(f"\n--- {panel_name} | {role} | N_total={total:,} ---")
    panel_set = set(panel_drugs)

    if role=='PS_SS':
        p = pairs[pairs['role_cod'].isin(['PS','SS']) & pairs['ingredient'].isin(panel_set)]
        d = drug_adol[drug_adol['role_cod'].isin(['PS','SS'])]
    else:
        p = pairs[pairs['ingredient'].isin(panel_set)]
        d = drug_adol.copy()

    # Case sets
    all_ids = set(adol['caseid'])
    drug_cases = d.groupby('ingredient')['caseid'].apply(set).to_dict()
    evt_cases  = reac_adol.groupby('pt')['caseid'].apply(set).to_dict()
    pair_cases = p.groupby(['ingredient','pt'])['caseid'].apply(set).to_dict()

    rows=[]
    unique = p[['ingredient','pt']].drop_duplicates()
    print(f"  Computing {len(unique):,} pairs...")
    for _,(ing,pt) in unique.iterrows():
        a_ = len(pair_cases.get((ing,pt), set()))
        if a_ < 3: continue
        d_ = drug_cases.get(ing, set())
        e_ = evt_cases.get(pt, set())
        b_ = len(d_ - e_)
        c_ = len(e_ - d_)
        dd = total - a_ - b_ - c_
        if dd < 0: continue

        rv,rl,ru = ror(a_,b_,c_,dd)
        pv,pl,pu = prr(a_,b_,c_,dd)
        iv,i0   = ic_bcpnn(a_,b_,c_,dd)
        ev,e0   = ebgm(a_,b_,c_,dd)

        rs = rl is not None and rl>1
        ps = pl is not None and pl>1
        is_ = i0 is not None and i0>0
        es = e0 is not None and e0>1
        nm = sum([rs,ps,is_,es])

        rows.append({
            'panel':panel_name,'drug':ing,'event':pt,'role_filter':role,
            'N':a_,'a':a_,'b':b_,'c':c_,'d':dd,
            'ROR':round(rv,4) if rv else None,
            'ROR_lower':round(rl,4) if rl else None,
            'ROR_upper':round(ru,4) if ru else None,
            'ROR_signal':rs,
            'PRR':round(pv,4) if pv else None,
            'PRR_lower':round(pl,4) if pl else None,
            'PRR_upper':round(pu,4) if pu else None,
            'PRR_signal':ps,
            'IC':round(iv,4) if iv else None,
            'IC025':round(i0,4) if i0 else None,
            'IC_signal':is_,
            'EBGM':round(ev,4) if ev else None,
            'EB05':round(e0,4) if e0 else None,
            'EBGM_signal':es,
            'n_metrics':nm,'four_metric':nm==4,
        })

    df=pd.DataFrame(rows)
    if len(df)>0:
        print(f"  Pairs N>=3: {len(df):,}  ROR signals: {df['ROR_signal'].sum()}  4-metric: {df['four_metric'].sum()}")
    else:
        print("  No pairs with N>=3")
    return df

def task_b(adol, drug_adol, reac_adol, pairs):
    print("\n"+"="*68)
    print("TASK B: DISPROPORTIONALITY")
    print("="*68)
    panels=[
        (OBESITY_BROAD,'Broad Obesity'),
        (DIABETES_BROAD,'Broad Diabetes'),
        (OBESITY_SEL,'Selected Obesity'),
        (DIABETES_SEL,'Selected Diabetes'),
    ]
    all_dfs=[]
    for drugs,name in panels:
        for role in ['PS_SS','ALL']:
            df=run_disprop(adol,drug_adol,reac_adol,pairs,drugs,name,role)
            all_dfs.append(df)
    combined=pd.concat(all_dfs,ignore_index=True) if all_dfs else pd.DataFrame()
    combined.to_csv(os.path.join(OUTPUT_DIR,'disproportionality_all.csv'),index=False)
    combined.to_excel(os.path.join(OUTPUT_DIR,'ESM_8_Complete_Signal_Detection.xlsx'),index=False)
    print(f"\nAll disproportionality results saved ({len(combined):,} rows)")
    return combined


# ===================================================================
# TABLES 1 & 4
# ===================================================================
def gen_table1(adol, drug_adol):
    print("\n"+"="*68+"  \nTABLE 1: Drug case counts")
    rows=[]
    for panel_drugs,panel_name in [(OBESITY_BROAD,'Broad Obesity'),
                                   (DIABETES_BROAD,'Broad Diabetes'),
                                   (OBESITY_SEL,'Selected Obesity'),
                                   (DIABETES_SEL,'Selected Diabetes')]:
        for d in panel_drugs:
            ps = drug_adol[(drug_adol['ingredient']==d)&(drug_adol['role_cod'].isin(['PS','SS']))]['caseid'].nunique()
            al = drug_adol[drug_adol['ingredient']==d]['caseid'].nunique()
            rows.append({'panel':panel_name,'drug':d,'cases_PS_SS':ps,'cases_all_roles':al})
            if ps>0: print(f"  {panel_name} | {d}: {ps} PS/SS, {al} all")
    df=pd.DataFrame(rows)
    df.to_csv(os.path.join(OUTPUT_DIR,'table1_drug_counts.csv'),index=False)
    zero=df[df['cases_PS_SS']==0]['drug'].unique()
    if len(zero): print(f"  Drugs with 0 PS/SS (remove from Table 1): {list(zero)}")
    return df

def gen_table4(adol, drug_adol, outc_adol):
    print("\n"+"="*68+"\nTABLE 4: Demographics per panel")
    rows=[]
    for panel_drugs,panel_name in [(OBESITY_BROAD,'Broad Obesity'),(DIABETES_BROAD,'Broad Diabetes')]:
        pset=set(panel_drugs)
        ids=set(drug_adol[(drug_adol['ingredient'].isin(pset))&(drug_adol['role_cod'].isin(['PS','SS']))]['caseid'])
        sub=adol[adol['caseid'].isin(ids)]
        n=len(sub)
        if n==0: continue
        drole=drug_adol[(drug_adol['caseid'].isin(ids))&(drug_adol['ingredient'].isin(pset))]
        rc=drole['role_cod'].value_counts()
        outc=outc_adol[outc_adol['caseid'].isin(ids)]['outc_cod'].value_counts()
        row={'panel':panel_name,'N':n,
             'age_mean':round(sub['age_years'].mean(),2),
             'age_sd':round(sub['age_years'].std(),2),
             'female_pct':round((sub['sex']=='F').mean()*100,1),
             'PS_pct':round(rc.get('PS',0)/len(drole)*100,1) if len(drole) else 0,
             'SS_pct':round(rc.get('SS',0)/len(drole)*100,1) if len(drole) else 0,
             'C_pct': round(rc.get('C',0)/len(drole)*100,1) if len(drole) else 0,
        }
        for code in ['DE','LT','HO','DS','CA','RI','OT']:
            row[f'outc_{code}']=int(outc.get(code,0))
        rows.append(row)
        print(f"  {panel_name}: N={n} | age {row['age_mean']}±{row['age_sd']} | PS {row['PS_pct']}% C {row['C_pct']}%")
    df=pd.DataFrame(rows)
    df.to_csv(os.path.join(OUTPUT_DIR,'table4_demographics.csv'),index=False)
    return df


# ===================================================================
# TASK D: ACCEPTANCE TESTS
# ===================================================================
def task_d(adol, pairs, disprop, total):
    print("\n"+"="*68+"\nTASK D: ACCEPTANCE TESTS")
    res={}

    ac=adol['age_years'].value_counts(normalize=True)
    pct=ac.iloc[0]*100
    res['T1_age_no_single_dominates_25pct']={'passed':pct<=25,'max_pct':round(pct,2),'at_age':float(ac.index[0])}
    print(f"  T1 age dominance: max {pct:.1f}% — {'PASS' if pct<=25 else 'FAIL'}")

    oor=((adol['age_years']<12)|(adol['age_years']>=18)).sum()
    res['T2_age_in_range']={'passed':int(oor)==0,'out_of_range':int(oor)}
    print(f"  T2 age range: {oor} out-of-range — {'PASS' if oor==0 else 'FAIL'}")

    if len(pairs)>0:
        nd=pairs.duplicated(subset=['caseid','ingredient','pt']).sum()
        res['T3_no_duplicate_triples']={'passed':int(nd)==0,'duplicates':int(nd)}
        print(f"  T3 duplicates: {nd} — {'PASS' if nd==0 else 'FAIL'}")

        if len(disprop)>0:
            sums=disprop['a']+disprop['b']+disprop['c']+disprop['d']
            ok=(sums==total).all()
            res['T4_2x2_sums_correct']={'passed':bool(ok),'expected':total}
            print(f"  T4 2×2 sums={total}: {'PASS' if ok else 'FAIL'}")

            sig=disprop[disprop['ROR_signal']]
            mn=int(sig['N'].min()) if len(sig)>0 else None
            res['T5_signal_N_ge_3']={'passed':mn is None or mn>=3,'min_N':mn}
            print(f"  T5 min signal N: {mn} — {'PASS' if (mn is None or mn>=3) else 'FAIL'}")

    n_pass=sum(1 for v in res.values() if v['passed'])
    print(f"\n  RESULT: {n_pass}/{len(res)} PASSED")
    with open(os.path.join(OUTPUT_DIR,'acceptance_tests.json'),'w') as f:
        json.dump(res,f,indent=2)
    return res


# ===================================================================
# MAIN
# ===================================================================
if __name__=='__main__':
    print("FAERS Adolescent Pipeline — Memory-Efficient Build")
    print(f"Output: {OUTPUT_DIR}\n")

    adol, drug_adol, reac_adol, outc_adol, pairs, stats = task_a()
    t1  = gen_table1(adol, drug_adol)
    t4  = gen_table4(adol, drug_adol, outc_adol)
    dis = task_b(adol, drug_adol, reac_adol, pairs)
    tst = task_d(adol, pairs, dis, stats['adolescent_cases'])

    # Master summary
    ps=dis[dis['role_filter']=='PS_SS'] if len(dis)>0 else pd.DataFrame()
    summary={'cohort':stats,'signals':{}}
    for pname in (ps['panel'].unique() if len(ps)>0 else []):
        sub=ps[ps['panel']==pname]
        summary['signals'][pname]={
            'ror_signals':int(sub['ROR_signal'].sum()),
            'four_metric':int(sub['four_metric'].sum()),
            'pairs_N_ge3':int(len(sub))
        }
    with open(os.path.join(OUTPUT_DIR,'pipeline_summary.json'),'w') as f:
        json.dump(summary,f,indent=2)

    print(f"\n{'='*68}")
    print("PIPELINE COMPLETE")
    print(f"  Adolescent cases : {stats['adolescent_cases']:,}")
    print(f"  Age : {stats['age_mean']} ± {stats['age_sd']}")
    for pn,sv in summary['signals'].items():
        print(f"  {pn} (PS/SS): {sv['ror_signals']} ROR signals, {sv['four_metric']} 4-metric")
    print(f"All outputs in: {OUTPUT_DIR}")
