"""Attempt to falsify the report's main descriptive conclusions.

Run: python analysis_review/stress_test_findings.py
Inputs and earlier outputs are read-only. Results go into outputs/sanity_checks.
Perturbations, permutations, and sequence components are sensitivity diagnostics,
not experimental error models or population confidence intervals.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import os
import platform
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd
import scipy
from scipy.stats import kendalltau
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/sanity_checks"
OUT.mkdir(parents=True, exist_ok=True)
MAP = {"APEX 1.0": "apex_1_0_E._coli_ATCC11775", "APEX 1.1": "apex_1_1_E._coli_ATCC_11775",
       "MBC-Attention": "mbc_attention_E._coli_unspecified_strain", "AMPredictor": "ampredictor_E._coli_unspecified_strain"}
PRIMARY = ["APEX 1.1", "MBC-Attention", "AMPredictor"]
CHECKS = []


def save(rows, name):
    pd.DataFrame(rows).to_csv(OUT / (name + ".csv"), index=False)


def check(name, condition):
    CHECKS.append({"check": name, "passed": bool(condition)})


def shortlist(z, k):
    return {m: set(z[m].nsmallest(k).index) for m in z}


def overlap_rows(z, name, k=None):
    if k is None:
        k = int(np.ceil(.01 * len(z)))
    ss = shortlist(z.sort_index(), k)
    return [{"population": name, "left": a, "right": b, "n": len(z), "k": k,
             "intersection": len(ss[a] & ss[b]), "shared_fraction": len(ss[a] & ss[b]) / k}
            for a, b in combinations(z, 2)]


def weighted_percentile(x, weight):
    order = np.argsort(x, kind="stable")
    xx, ww = x[order], weight[order]
    # Midrank weighted ECDF, sharing identical values exactly.
    unique, first, counts = np.unique(xx, return_index=True, return_counts=True)
    totals = np.add.reduceat(ww, first)
    middle = (np.cumsum(totals) - totals / 2) / totals.sum()
    result = np.empty(len(x))
    result[order] = np.repeat(middle, counts)
    return result


def main():
    started = time.monotonic()
    path = ROOT / "df_damp_195625x57.tsv"
    raw = pd.read_csv(path, sep="\t").set_index("sequence", verify_integrity=True)
    source = pd.read_csv(ROOT / "scratchpad/source_sequences.tsv", sep="\t").set_index("sequence", verify_integrity=True).source
    v = raw[list(MAP.values())].rename(columns={c:m for m,c in MAP.items()}).dropna().sort_index()
    z = np.log2(v)
    f = pd.DataFrame({"source": source.reindex(z.index), "length": z.index.str.len()}, index=z.index)
    f["has_cysteine"] = z.index.str.contains("C")
    n = len(z); k = int(np.ceil(.01*n)); ss = shortlist(z, k)
    primary = z[PRIMARY]
    a = primary.to_numpy()
    ranks = primary.rank(pct=True).to_numpy()
    spread = np.ptp(a, axis=1)
    lo, hi = np.quantile(spread, [.05,.95])
    low = spread <= lo; high = spread >= hi
    group = np.select([low, high], ["low", "high"], default="middle")
    print("Checking numeric invariants and cutoffs", flush=True)
    check("same complete cases", n == 194947)
    check("positive finite complete predictions", np.isfinite(v).all().all() and v.gt(0).all().all())
    check("log base cannot change raw-tail membership", np.array_equal(np.log10(v[PRIMARY]).max(axis=1).sub(np.log10(v[PRIMARY]).min(axis=1)).le(lo / np.log2(10)).to_numpy(), low))
    check("common concentration-unit change leaves spreads unchanged", np.allclose(np.ptp(np.log2(v[PRIMARY]*1e-6),axis=1), spread))
    centered = primary - primary.median()
    check("centering preserves all primary shortlists", all(set(centered[m].nsmallest(k).index) == ss[m] for m in PRIMARY))
    shuffled = z.sample(frac=1, random_state=704)
    check("row order cannot change primary shortlists", all(set(shuffled[m].nsmallest(k).index)==ss[m] for m in PRIMARY))
    tie_rows = []
    for model in z:
        zz = z[model].to_numpy(); cutoff = np.partition(zz, k-1)[k-1]
        independent = set(z.index[np.argpartition(zz,k-1)[:k]])
        check("numpy independently reproduces shortlist " + model, independent == ss[model])
        tie_rows.append({"model": model, "k": k, "cutoff_uM": float(2**cutoff),
                         "n_strictly_below_cutoff": int((zz<cutoff).sum()), "n_equal_cutoff": int((zz==cutoff).sum()),
                         "n_within_0_1_log2_of_cutoff": int((abs(zz-cutoff)<=.1).sum()),
                         "n_within_0_01_log2_of_cutoff": int((abs(zz-cutoff)<=.01).sum())})
    save(tie_rows,"shortlist_ties_and_boundaries")
    # Directly enumerate a sample's pair orders, including deliberately created ties.
    r = np.random.default_rng(704)
    tiny = z.iloc[r.choice(n, size=250, replace=False)].round(1)
    for left,right in combinations(tiny,2):
        x,y=tiny[left].to_numpy(),tiny[right].to_numpy()
        dx,dy=np.subtract.outer(x,x),np.subtract.outer(y,y)
        upper=np.triu(np.ones(dx.shape,dtype=bool),k=1)
        comparable=upper & (dx!=0) & (dy!=0)
        direct=float((dx[comparable]*dy[comparable]<0).mean())
        total=len(x)*(len(x)-1)/2
        tx=float((upper & (dx==0)).sum());ty=float((upper & (dy==0)).sum())
        tau=kendalltau(x,y).statistic
        inferred=(comparable.sum()-tau*np.sqrt((total-tx)*(total-ty)))/(2*comparable.sum())
        check("Kendall reversal formula matches explicit pairs "+left+"/"+right, np.isclose(direct,inferred))

    # Hypothetical bounded perturbations; magnitude is not measured model noise.
    print("Testing shortlist perturbations and substantial pair reversals",flush=True)
    jitter=[]
    rng=np.random.default_rng(705)
    for bound in [.01,.1,.25,.5,1.]:
        for repeat in range(25):
            noise=rng.uniform(-bound,bound,size=z.shape)
            selections=[]
            for j,model in enumerate(z):
                idx=np.argpartition(z.iloc[:,j].to_numpy()+noise[:,j],k-1)[:k]
                selection=set(z.index[idx]);selections.append(selection)
                jitter.append({"kind":"same_model_original_vs_perturbed","bound_log2":bound,"repeat":repeat,
                               "left":model,"right":model,"shared_fraction":len(selection&ss[model])/k})
            for i,j in combinations(range(4),2):
                jitter.append({"kind":"two_perturbed_models","bound_log2":bound,"repeat":repeat,
                               "left":z.columns[i],"right":z.columns[j],"shared_fraction":len(selections[i]&selections[j])/k})
    jitter=pd.DataFrame(jitter);save(jitter,"bounded_score_perturbations")
    save(jitter.groupby(["kind","bound_log2","left","right"]).shared_fraction.agg(["min","median","max"]).reset_index(),"bounded_score_perturbation_summary")
    directed=[]
    rank_all=z.rank(pct=True)
    for left,right in [(l,r) for l in z for r in z if l!=r]:
        idx=sorted(ss[left]); cutoff_left=z[left].nsmallest(k).max();cutoff_right=z[right].nsmallest(k).max()
        directed.append({"selected_by":left,"assessed_by":right,"k":k,
                         "median_assessed_percentile":rank_all.loc[idx,right].median(),
                         "fraction_outside_assessed_best_10pct":rank_all.loc[idx,right].gt(.1).mean(),
                         "n_left_at_least_twofold_inside_and_right_twofold_outside":int(((z.loc[idx,left]<=cutoff_left-1)&(z.loc[idx,right]>=cutoff_right+1)).sum()),
                         "n_left_at_least_fourfold_inside_and_right_fourfold_outside":int(((z.loc[idx,left]<=cutoff_left-2)&(z.loc[idx,right]>=cutoff_right+2)).sum())})
    save(directed,"cross_model_shortlist_positions")
    # Sample unordered pairs uniformly with replacement, retaining both directions
    # in the sampling distribution. Report conditional margins and their frequency.
    rng=np.random.default_rng(706); pair_n=2_000_000
    ii=rng.integers(n,size=pair_n);jj=rng.integers(n,size=pair_n)
    keep=ii!=jj;ii,jj=ii[keep],jj[keep]
    margins=[]
    for left,right in combinations(z,2):
        dl=z[left].to_numpy()[ii]-z[left].to_numpy()[jj]
        dr=z[right].to_numpy()[ii]-z[right].to_numpy()[jj]
        for margin in [0,.5,1.,2.]:
            ok=(abs(dl)>margin)&(abs(dr)>margin)
            margins.append({"left":left,"right":right,"required_log2_difference_each_model":margin,
                            "sampled_pairs":len(ii),"n_comparable":int(ok.sum()),"eligible_fraction":ok.mean(),
                            "opposite_fraction_among_eligible":(dl[ok]*dr[ok]<0).mean()})
    save(margins,"order_reversals_with_effect_margins")

    print("Testing sources, missingness, and budgets",flush=True)
    pop_rows=[]
    populations={"all":np.ones(n,dtype=bool),
                 "omit_Archaeomics_and_Extinctomics":~f.source.isin(["Archaeomics","Extinctomics"]).to_numpy(),
                 "omit_database_collection":~f.source.eq("DRAMP_APD_DBAASP").to_numpy(),
                 "cysteine_free_8_to_39":(f.length.between(8,39)&~f.has_cysteine).to_numpy()}
    training=set(pd.read_csv(ROOT/"scratchpad/replay/APEXDB.csv").Sequence)
    for line in (ROOT/"scratchpad/replay/ampredictor_train.txt").read_text().splitlines():
        bits=line.replace("\\t","\t").split("\t")
        if len(bits)>=2:training.add(bits[1].strip())
    populations["omit_exact_cached_training_matches"]=~z.index.isin(training)
    for src in sorted(f.source.unique()):
        populations["leave_out_"+src]=~f.source.eq(src).to_numpy()
    for name,mask in populations.items():
        pop_rows.extend(overlap_rows(z.loc[mask],name))
    save(pop_rows,"population_shortlist_sensitivity")
    source_budget={m:set() for m in z}
    for src,ix in f.groupby("source").groups.items():
        b=z.loc[ix];kk=int(np.ceil(.01*len(b)))
        for m,sel in shortlist(b,kk).items():source_budget[m].update(sel)
    save([{"left":l,"right":r,"budget":len(source_budget[l]),"intersection":len(source_budget[l]&source_budget[r]),
           "shared_fraction":len(source_budget[l]&source_budget[r])/len(source_budget[l])}
          for l,r in combinations(z,2)],"source_proportional_shortlists")
    budget_rows=[]
    for budget in [100,250,500,1000,1950,5000,10000,20000]:
        budget_rows.extend(overlap_rows(z,"all",budget))
    save(budget_rows,"absolute_screening_budgets")
    # Conservative bounds when adding rows absent from the pairwise universe.
    # Each model can lose old shortlist entries only to unknown own scores or
    # newly admitted known scores below its old cutoff. Old-old intersection
    # cannot increase; added rows can create at most min(replacements) new matches.
    miss=[];completion_rng=np.random.default_rng(708)
    for left,right in combinations(z,2):
        union_valid=raw[[MAP[left],MAP[right]]].notna().all(axis=1)
        missing_l=int(raw[MAP[left]].isna().sum());missing_r=int(raw[MAP[right]].isna().sum())
        ov=len(ss[left]&ss[right])
        # Additional universe rows can include rows whose other (third) model is
        # missing. Pairwise complete-case comparison below uses all available rows.
        vv=raw.loc[union_valid,[MAP[left],MAP[right]]].rename(columns={MAP[left]:left,MAP[right]:right})
        actual=overlap_rows(np.log2(vv),"pairwise_available",k)[0]
        cutoff_l=vv[left].nsmallest(k).max();cutoff_r=vv[right].nsmallest(k).max()
        known_new_l=int(((~union_valid)&raw[MAP[left]].notna()&raw[MAP[left]].le(cutoff_l)).sum())
        known_new_r=int(((~union_valid)&raw[MAP[right]].notna()&raw[MAP[right]].le(cutoff_r)).sum())
        replace_l=min(k,missing_l+known_new_l);replace_r=min(k,missing_r+known_new_r)
        lower=max(0,actual["intersection"]-replace_l-replace_r)
        upper=min(k,actual["intersection"]+min(replace_l,replace_r))
        completion_counts=[]
        full_logs=np.log2(raw[[MAP[left],MAP[right]]].to_numpy())
        holes=np.isnan(full_logs)
        for trial in range(50):
            completed=full_logs.copy()
            completed[holes]=completion_rng.uniform(-20,20,size=int(holes.sum()))
            il=np.argpartition(completed[:,0],k-1)[:k];ir=np.argpartition(completed[:,1],k-1)[:k]
            completion_counts.append(len(set(il)&set(ir)))
        check("missing-score completions respect derived bounds "+left+"/"+right,all(lower<=count<=upper for count in completion_counts))
        miss.append({"left":left,"right":right,"missing_left":missing_l,"missing_right":missing_r,
                     "pairwise_available_n":len(vv),"common_universe_intersection":ov,
                     "pairwise_available_intersection_same_k":actual["intersection"],
                     "max_possible_replacements_left":replace_l,"max_possible_replacements_right":replace_r,
                     "adversarial_full_universe_min_intersection":lower,
                     "adversarial_full_universe_max_intersection":upper,
                     "smallest_intersection_in_50_random_completions":min(completion_counts),
                     "largest_intersection_in_50_random_completions":max(completion_counts),
                     "bound_scope":"fixed k=1950; arbitrary finite missing scores; loose count bound, not probability"})
    save(miss,"missing_prediction_bounds")

    print("Testing whether agreement-tail geometry needs peptide-specific alignment",flush=True)
    # Re-rank within source or source/length, and use equally weighted sources.
    local_rows=[]
    source_weight=(1/f.groupby("source").source.transform("size")).to_numpy()
    rank_defs={"global":ranks,
               "within_source":primary.groupby(f.source).rank(pct=True).to_numpy(),
               "within_source_exact_length":primary.groupby([f.source,f.length]).rank(pct=True).to_numpy(),
               "equal_source_weight_ECDF":np.column_stack([weighted_percentile(a[:,j],source_weight) for j in range(3)])}
    for name,rr in rank_defs.items():
        for g in ["low","middle","high"]:
            mask=group==g
            local_rows.append({"rank_reference":name,"raw_group":g,"n":int(mask.sum()),
                               "median_rank_spread":np.median(np.ptp(rr[mask],axis=1))})
    save(local_rows,"rank_reference_sensitivity")
    # Permutations preserve marginal model outputs; source-length permutations
    # additionally preserve their observed association with source and exact length.
    # They destroy alignment of predictions to the same peptide within each stratum.
    permutations=[];rng=np.random.default_rng(707)
    strata=list(f.groupby(["source","length"]).indices.values())
    def summarize_arrays(b,rr,condition,repeat):
        sp=np.ptp(b,axis=1);qlo,qhi=np.quantile(sp,[.05,.95]);ml=sp<=qlo;mh=sp>=qhi
        flags=rr <= k/n
        row={"condition":condition,"repeat":repeat,"median_raw_fold_spread":2**np.median(sp),
             "raw_low_median_rank_spread":np.median(np.ptp(rr[ml],axis=1)),
             "raw_high_median_rank_spread":np.median(np.ptp(rr[mh],axis=1)),
             "high_minus_low_median_length":np.median(f.length.to_numpy()[mh])-np.median(f.length.to_numpy()[ml])}
        for i,j in combinations(range(3),2):row[PRIMARY[i]+" / "+PRIMARY[j]+" shortlist_count"]=int((flags[:,i]&flags[:,j]).sum())
        return row
    permutations.append(summarize_arrays(a,ranks,"observed",-1))
    for condition in ["independent_global_permutations","independent_source_length_permutations"]:
        for repeat in range(30):
            orders=[]
            for j in range(3):
                if condition=="independent_global_permutations":idx=rng.permutation(n)
                else:
                    idx=np.arange(n)
                    for members in strata:idx[members]=rng.permutation(members)
                orders.append(idx)
            b=np.column_stack([a[orders[j],j] for j in range(3)])
            rr=np.column_stack([ranks[orders[j],j] for j in range(3)])
            permutations.append(summarize_arrays(b,rr,condition,repeat))
    perm=pd.DataFrame(permutations);save(perm,"agreement_geometry_permutations")
    perm_summary=perm.groupby("condition").agg({c:["min","median","max"] for c in perm.select_dtypes("number") if c!="repeat"})
    perm_summary.columns=[name+"_"+stat for name,stat in perm_summary.columns]
    save(perm_summary.reset_index(),"agreement_geometry_permutation_summary")

    print("Constructing broader sequence components using exact signatures",flush=True)
    # Work one length and deletion position at a time. The equal-position deletion
    # signature identifies exact Hamming-distance-one substitutions, without an
    # all-pairs alignment or a large global signature dictionary.
    sequences=list(raw.index);lookup={s:i for i,s in enumerate(sequences)}
    parent=np.arange(len(sequences));size=np.ones(len(sequences),dtype=int)
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(i,j):
        ri,rj=root(i),root(j)
        if ri==rj:return
        if size[ri]<size[rj]:ri,rj=rj,ri
        parent[rj]=ri;size[ri]+=size[rj]
    terminal_edges=0;substitution_edges=0
    length_groups={}
    for i,s in enumerate(sequences):
        length_groups.setdefault(len(s),[]).append(i)
        if len(s)>=9:
            for t in {s[1:],s[:-1]}:
                if t in lookup:union(i,lookup[t]);terminal_edges+=1
    for length,indices in sorted(length_groups.items()):
        if length<8:continue
        for pos in range(length):
            seen={}
            for i in indices:
                s=sequences[i];signature=s[:pos]+s[pos+1:]
                if signature in seen:union(i,seen[signature]);substitution_edges+=1
                else:seen[signature]=i
    check("terminal edges reproduce previous count",terminal_edges==92233)
    component=pd.Series({s:root(i) for i,s in enumerate(sequences)})
    comp=component.reindex(z.index)
    counts=comp.value_counts()
    sizes=counts.reindex(comp).to_numpy()
    save(pd.DataFrame({"sequence":z.index,"component":comp.to_numpy(),"complete_case_component_size":sizes}),"broader_sequence_components")
    component_rows=[]
    for left,right in combinations(z,2):
        ca=set(comp.loc[sorted(ss[left])]);cb=set(comp.loc[sorted(ss[right])]);inter=len(ca&cb)
        component_rows.append({"left":left,"right":right,"unique_selected_components_left":len(ca),
                               "unique_selected_components_right":len(cb),"shared_components":inter,
                               "component_recall_left":inter/len(ca),"component_recall_right":inter/len(cb),
                               "component_jaccard":inter/len(ca|cb)})
    save(component_rows,"original_shortlists_component_overlap")
    reps_rows=[]
    for seed in range(12):
        hashes=np.array([hashlib.sha256((str(seed)+":"+s).encode()).hexdigest() for s in z.index])
        order=np.argsort(hashes)
        c_sorted=comp.to_numpy()[order];_,first=np.unique(c_sorted,return_index=True)
        representatives=z.iloc[order[first]]
        for record in overlap_rows(representatives,"broader_component_representatives"):
            record["representative_seed"]=seed;reps_rows.append(record)
    reps=pd.DataFrame(reps_rows);save(reps,"broader_component_representative_shortlists")
    save(reps.groupby(["left","right"]).shared_fraction.agg(["min","median","max"]).reset_index(),"broader_component_representative_summary")
    family_medians=z.groupby(comp.to_numpy()).median()
    save(overlap_rows(family_medians,"component_median_prediction"),"component_median_shortlists")
    save([{"n_full_input_components":int(component.nunique()),"n_complete_case_components":int(comp.nunique()),
           "n_nontrivial_complete_case_components":int(counts.gt(1).sum()),"largest_component_complete_cases":int(counts.max()),
           "fraction_complete_cases_in_nontrivial_components":float((sizes>1).mean()),
           "terminal_edges":terminal_edges,"substitution_union_edges":substitution_edges,
           "scope":"single-linkage terminal extension plus equal-length one-substitution edges; minimum shorter/substitution length 8"}],"broader_component_summary")

    print("Auditing provisional experimental labels with interval and censoring losses",flush=True)
    e=pd.read_csv(ROOT/"analysis_review/outputs/archaeasin_provisional_ecoli_labels.csv")
    numeric=e.MIC_uM.notna();audit_labels=[];loo=[]
    for model in MAP:
        p=np.log2(e[model].to_numpy());observed=np.log2(e.MIC_uM.to_numpy())
        point_error=abs(p[numeric]-observed[numeric])
        interval_loss=np.maximum.reduce([observed[numeric]-1-p[numeric],p[numeric]-observed[numeric],np.zeros(numeric.sum())])
        inactive_lower_loss=np.maximum(6-p[~numeric],0)
        pseudo=np.where(numeric,observed,9.)
        audit_labels.append({"model":model,"n_numeric":int(numeric.sum()),"n_source_coded_inactive":int((~numeric).sum()),
                             "numeric_point_MAE_log2":point_error.mean(),"numeric_doubling_interval_distance_log2":interval_loss.mean(),
                             "all79_distance_to_intervals_conditional_inactive_gt64_log2":np.r_[interval_loss,inactive_lower_loss].mean(),
                             "all79_MAE_if_wrongly_treat_512_as_exact":abs(p-pseudo).mean(),
                             "n_numeric_predictions_below_interval":int((p[numeric]<observed[numeric]-1).sum()),
                             "n_numeric_predictions_inside_interval":int(((p[numeric]>=observed[numeric]-1)&(p[numeric]<=observed[numeric])).sum()),
                             "n_numeric_predictions_above_interval":int((p[numeric]>observed[numeric]).sum())})
        for j in range(len(point_error)):
            loo.append({"model":model,"omitted_numeric_position":j,"remaining_numeric_MAE_log2":np.delete(point_error,j).mean()})
    save(audit_labels,"experimental_label_loss_sensitivity");save(loo,"experimental_leave_one_numeric_out")
    manifest=json.loads((ROOT/"analysis_review/outputs/archaeasin_curation_manifest.json").read_text())
    for file in manifest["files"]:
        check("downloaded source checksum "+file["name"],hashlib.sha256((ROOT/"analysis_review/inputs"/file["name"]).read_bytes()).hexdigest()==file["sha256"])
    idmap=pd.read_csv(ROOT/"analysis_review/outputs/archaeasin_provisional_id_map.csv")
    check("inferred numbered identities remain unique",idmap.sequence.is_unique)
    check("all fingerprint leave-target-out checks passed",bool(idmap.all_11_leave_target_out_matches_agree.all()))
    check("one chemistry discrepancy excluded",len(e)==79 and int(idmap.exact_sequence_in_synthesis_certificates.sum())==79)

    # Fixed sensitivity based on primary published label coding, not an offset
    # estimated to improve agreement on this cohort. It is NOT a calibration.
    # Original APEX methods: inactive MIC >128 coded 140 uM.
    # APEX 1.1 methods: inactive above tested maximum coded 512 uM.
    print("Testing the published APEX inactive-label scale change",flush=True)
    metadata=pd.DataFrame(json.loads((ROOT/"scratchpad/damp_metadata.json").read_text())).set_index("column")
    older=metadata[metadata.model.eq("apex_1_0")]
    newer=metadata[metadata.model.eq("apex_1_1")]
    coding_factor=512/140;coding_rows=[]
    for new_col,row in newer.iterrows():
        if bool(row.alias_unverified):continue
        old_col=older.index[older.strain_index.eq(row.strain_index)][0]
        if bool(older.loc[old_col,"alias_unverified"]):continue
        ratio=raw[new_col]/raw[old_col];residual=np.log2(ratio/coding_factor)
        coding_rows.append({"target":row.target,"n":len(ratio),"median_raw_v11_v10_ratio":ratio.median(),
                            "fixed_pseudolabel_factor":coding_factor,"median_ratio_after_fixed_factor":(ratio/coding_factor).median(),
                            "raw_within_twofold":np.log2(ratio).abs().le(1).mean(),
                            "within_twofold_after_fixed_factor":residual.abs().le(1).mean(),
                            "median_absolute_log2_difference_after_fixed_factor":residual.abs().median()})
    save(coding_rows,"apex_pseudolabel_scale_sensitivity")
    check("eight unambiguous targets in label-coding sensitivity",len(coding_rows)==8)
    xx=raw[MAP["APEX 1.0"]];yy=raw[MAP["APEX 1.1"]];ratio=yy/xx;coding_bins=[]
    for lower,upper in [(0,8),(8,16),(16,32),(32,64),(64,128),(128,np.inf)]:
        mask=xx.gt(lower)&xx.le(upper)
        coding_bins.append({"APEX_v10_prediction_lower_uM":lower,"APEX_v10_prediction_upper_uM":upper,
                            "n":int(mask.sum()),"median_v10_prediction_uM":xx[mask].median(),
                            "median_v11_v10_ratio":ratio[mask].median(),
                            "median_ratio_after_fixed_pseudolabel_factor":(ratio[mask]/coding_factor).median(),
                            "within_twofold_after_fixed_factor":np.log2(ratio[mask]/coding_factor).abs().le(1).mean()})
    save(coding_bins,"apex_pseudolabel_scale_by_prediction_regime")
    (OUT/"apex_label_coding_sources.json").write_text(json.dumps({
        "original_APEX_methods":"https://pmc.ncbi.nlm.nih.gov/articles/PMC11310081/",
        "original_APEX_article":"https://doi.org/10.1038/s41551-024-01201-x",
        "original_inactive_threshold_uM":128,"original_inactive_pseudolabel_uM":140,
        "APEX_1_1_methods":"https://www.nature.com/articles/s41564-025-02061-0",
        "APEX_1_1_inactive_pseudolabel_uM":512,"ratio":coding_factor,
        "interpretation":"Plausible scale contributor; training corpora and tasks also changed. No causal attribution or calibrated conversion is established."},indent=2)+"\n")
    save(CHECKS,"sanity_check_results")

    # Two concise figures documenting checks that most directly affect interpretation.
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    for model,color in zip(z,["#9c7770","#225b78","#c16936","#287a65"]):
        subset=jitter[(jitter.kind=="same_model_original_vs_perturbed")&jitter.left.eq(model)]
        b=subset.groupby("bound_log2").shared_fraction.agg(["min","median","max"])
        axes[0].plot(b.index,100*b["median"],marker="o",label=model,color=color)
        axes[0].fill_between(b.index,100*b["min"],100*b["max"],color=color,alpha=.12)
    axes[0].set(xlabel="Assumed independent uniform perturbation bound (log2)",ylabel="Original 1% shortlist retained (%)",ylim=(0,102),title="A. Same-model stability under hypothetical perturbations")
    axes[0].legend(fontsize=8)
    rows=pd.DataFrame(budget_rows)
    for left,right,color in [("APEX 1.1","MBC-Attention","#225b78"),("APEX 1.1","AMPredictor","#287a65"),("MBC-Attention","AMPredictor","#c16936")]:
        b=rows[rows.left.eq(left)&rows.right.eq(right)]
        axes[1].plot(b.k,100*b.shared_fraction,marker="o",label=left+" / "+right,color=color)
    axes[1].set(xscale="log",xlabel="Peptides selected by each model",ylabel="Shared list fraction (%)",title="B. Overlap depends on the screening budget")
    axes[1].legend(fontsize=7.5)
    fig.suptitle("Small bounded score perturbations do not close the cross-model shortlist gap",fontsize=13,weight="bold")
    for ext in ["png","pdf"]:fig.savefig(OUT/("figure_4_screening_stress_tests."+ext),dpi=180,bbox_inches="tight")
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    conditions=["independent_global_permutations","independent_source_length_permutations"]
    for i,condition in enumerate(conditions):
        b=perm[perm.condition.eq(condition)]
        axes[i].hist(100*b.raw_low_median_rank_spread,bins=10,color="#a8bac5",label="30 independent permutations")
        axes[i].axvline(100*np.median(np.ptp(ranks[low],axis=1)),color="#b04a41",lw=2,label="Observed")
        axes[i].set(xlabel="Raw-low group's median rank spread (percentage points)",ylabel="Permutation runs",title=["A. Preserve only model marginal distributions","B. Also preserve source and exact length"][i])
        axes[i].legend(fontsize=8)
    fig.suptitle("A large rank spread in the raw-low tail can arise from selection geometry",fontsize=13,weight="bold")
    for ext in ["png","pdf"]:fig.savefig(OUT/("figure_5_agreement_selection_geometry."+ext),dpi=180,bbox_inches="tight")
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.2,5),layout="constrained")
    bins=pd.DataFrame(coding_bins)
    ax.hexbin(np.log2(xx),np.log2(ratio),gridsize=70,bins="log",cmap="Blues",mincnt=1)
    ax.plot(np.log2(bins.median_v10_prediction_uM),np.log2(bins.median_v11_v10_ratio),"o-",color="#b04a41",label="Median within v1.0 prediction bands")
    ax.axhline(np.log2(coding_factor),ls="--",color="#287a65",label="Published pseudo-label ratio: 512 / 140 = 3.66")
    ax.set(xticks=np.log2([2,4,8,16,32,64,128]),xticklabels=[2,4,8,16,32,64,128],
           yticks=np.log2([.25,.5,1,2,4,8,16]),yticklabels=[.25,.5,1,2,4,8,16],
           xlabel="APEX v1.0 E. coli predicted MIC (µM; log scale)",
           ylabel="APEX v1.1 / v1.0 prediction ratio (fold; log scale)",
           title="A published label-coding change closely matches the global scale shift")
    ax.legend(fontsize=8,loc="lower left")
    for ext in ["png","pdf"]:fig.savefig(OUT/("figure_6_apex_label_coding_sensitivity."+ext),dpi=180,bbox_inches="tight")
    plt.close(fig)
    summary={"input_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"python":platform.python_version(),
             "numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,
             "n":n,"k":k,"seed_start":704,"permutation_repeats_each":30,"perturbation_repeats_each":25,
             "pair_sampling_n":len(ii),"family_representative_seeds":12,"checks":len(CHECKS),
             "passed":sum(r["passed"] for r in CHECKS),"failed":[r for r in CHECKS if not r["passed"]],
             "elapsed_seconds":time.monotonic()-started,
             "scope":"Exploratory falsification and sensitivity checks; no experimental uncertainty or population significance inferred."}
    (OUT/"stress_test_manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    assert all(r["passed"] for r in CHECKS),"See recorded sanity-check failures"


if __name__=="__main__":main()
