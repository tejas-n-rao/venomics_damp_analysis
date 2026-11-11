"""Absolute MIC agreement with divergent relative ranking, and sequence profiles.

Run: python analysis_review/similar_mic_different_ranks.py
Primary definition: all three exported MICs within twofold AND max-minus-min
within-model percentile >= 0.50, on the same 194,947 complete-case sequences.
Thresholds are exploratory descriptive definitions, not experimental tolerances.
"""
from pathlib import Path
from itertools import combinations, permutations
from collections import Counter
import hashlib
import json
import os
import platform
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"venomics-mpl"))
os.environ["MPLBACKEND"]="Agg"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"analysis_review/outputs/similar_mic_different_ranks"
OUT.mkdir(parents=True,exist_ok=True)
MAP={"APEX 1.1":"apex_1_1_E._coli_ATCC_11775","MBC-Attention":"mbc_attention_E._coli_unspecified_strain",
     "AMPredictor":"ampredictor_E._coli_unspecified_strain"}
CHECKS=[]


def save(frame,name):
    pd.DataFrame(frame).to_csv(OUT/(name+".csv"),index=False)


def check(name,condition):
    CHECKS.append({"check":name,"passed":bool(condition)})


def reference_percentile(sorted_values,x):
    left=np.searchsorted(sorted_values,x,side="left");right=np.searchsorted(sorted_values,x,side="right")
    # For observed values, exactly reproduce pandas' average tied ranks / n.
    return np.where(right>left,(left+1+right)/2,left)/len(sorted_values)


def main():
    start=time.monotonic();input_path=ROOT/"df_damp_195625x57.tsv"
    values=pd.read_csv(input_path,sep="\t",usecols=["sequence"]+list(MAP.values())).set_index("sequence",verify_integrity=True)
    values=values.rename(columns={v:k for k,v in MAP.items()}).dropna().sort_index()
    z=np.log2(values);ranks=values.rank(pct=True);n=len(values)
    spread=np.ptp(z,axis=1);rank_gap=np.ptp(ranks,axis=1)
    close=spread.le(1);target=close&rank_gap.ge(.5)
    source=pd.read_csv(ROOT/"scratchpad/source_sequences.tsv",sep="\t").set_index("sequence",verify_integrity=True).source.reindex(values.index)
    component=pd.read_csv(ROOT/"analysis_review/outputs/sanity_checks/broader_sequence_components.csv").set_index("sequence",verify_integrity=True).component.reindex(values.index)
    seq=values.index.to_series();length=seq.str.len()
    base=pd.DataFrame({"source":source,"length":length,"component":component,"median_log2_MIC":z.median(axis=1)},index=values.index)
    f=pd.DataFrame(index=values.index)
    f["length"]=length;f["basic_fraction"]=seq.str.count("[KR]")/length;f["acidic_fraction"]=seq.str.count("[DE]")/length
    f["charge_proxy"]=f.basic_fraction-f.acidic_fraction
    f["hydrophobic_fraction"]=seq.str.count("[AILMFWVY]")/length
    f["cysteine_fraction"]=seq.str.count("C")/length
    f["contains_cysteine"]=seq.str.contains("C").astype(float)
    f["no_D_or_E"]=~seq.str.contains("[DE]")
    f["positive_charge_proxy"]=f.charge_proxy.gt(0)
    f["entropy_bits"]=[-sum((k/len(s))*np.log2(k/len(s)) for k in Counter(s).values()) for s in seq]
    feature_names=["length","basic_fraction","acidic_fraction","charge_proxy","hydrophobic_fraction","cysteine_fraction",
                   "contains_cysteine","no_D_or_E","positive_charge_proxy","entropy_bits"]
    basic=seq.str.count("[KR]");hydrophobic=seq.str.count("[AILMFWVY]")
    basic_adj=np.array([sum(a in "KR" and b in "KR" for a,b in zip(s,s[1:])) for s in seq])
    hydro_adj=np.array([sum(a in "AILMFWVY" and b in "AILMFWVY" for a,b in zip(s,s[1:])) for s in seq])
    f["basic_adjacency_excess"]=(basic_adj-basic*(basic-1)/length)/(length-1)
    f["hydrophobic_adjacency_excess"]=(hydro_adj-hydrophobic*(hydrophobic-1)/length)/(length-1)
    for aa in "ACDEFGHIKLMNPQRSTVWY":f[aa+"_fraction"]=seq.str.count(aa)/length
    check("complete case denominator",n==194947)
    check("positive finite exports and source/component coverage",np.isfinite(z).all().all() and base.notna().all().all())
    check("direct ratio and log range agree",np.array_equal(values.max(axis=1).div(values.min(axis=1)).le(2),close))
    check("all subgroup cases are MIC-agreed",bool((~target|close).all()))
    for model in values:
        check("reference percentiles reproduce tied ranks "+model,
              np.allclose(reference_percentile(np.sort(values[model]),values[model]),ranks[model]))
    toy="KKALA";perm=list(permutations(toy))
    expected=2*(2-1)/len(toy)
    observed=np.mean([sum(a in "KR" and b in "KR" for a,b in zip(s,s[1:])) for s in perm])
    check("composition-preserving adjacency expectation",np.isclose(expected,observed))

    print("Quantifying MIC agreement and joint rank discordance",flush=True)
    prevalence=[]
    for tol in [1.1,1.25,1.5,2,4,8,16]:
        agreed=spread.le(np.log2(tol))
        for gap in [.25,.5,.75]:
            joint=agreed&rank_gap.ge(gap)
            prevalence.append({"maximum_fold_MIC_spread":tol,"minimum_rank_gap_percentage_points":100*gap,"n_complete":n,
                               "MIC_agreement_n":int(agreed.sum()),"MIC_agreement_fraction_of_all":agreed.mean(),
                               "joint_n":int(joint.sum()),"joint_fraction_of_all":joint.mean(),
                               "joint_fraction_of_MIC_agreement":joint.sum()/agreed.sum() if agreed.sum() else np.nan,
                               "median_rank_gap_among_MIC_agreement":rank_gap.loc[agreed].median()})
    save(prevalence,"prevalence_threshold_sensitivity")
    pairwise=[]
    for left,right in combinations(values,2):
        for tol in [2,4]:
            agreed=z[left].sub(z[right]).abs().le(np.log2(tol))
            gap=ranks[left].sub(ranks[right]).abs()
            pairwise.append({"left":left,"right":right,"maximum_fold_MIC_spread":tol,"n_common":n,
                             "MIC_agreement_n":int(agreed.sum()),"MIC_agreement_fraction":agreed.mean(),
                             "agreement_and_rank_gap_ge_50pp_n":int((agreed&gap.ge(.5)).sum()),
                             "rank_gap_ge_50pp_fraction_of_agreed":gap.loc[agreed].ge(.5).mean()})
    save(pairwise,"pairwise_agreement_and_ranks")
    percentiles=[]
    for mic in [8,16,32,64,128,256]:
        percentiles.append({"common_MIC_uM":mic,**{model+"_percentile":100*float(reference_percentile(np.sort(values[model]),mic)) for model in values}})
    save(percentiles,"same_value_different_reference_percentiles")
    profile=[]
    groups={"similar_MIC_large_rank_gap":target,"similar_MIC_smaller_rank_gap":close&~target,
            "all_similar_MIC":close,"rest_of_collection":~target,"all_complete":pd.Series(True,index=values.index)}
    for name,mask in groups.items():
        block=f.loc[mask].astype(float)
        for feature in f:
            profile.append({"group":name,"feature":feature,"n":len(block),"mean":block[feature].mean(),"median":block[feature].median(),
                            "q25":block[feature].quantile(.25),"q75":block[feature].quantile(.75)})
    save(profile,"sequence_profiles")
    by_source=[]
    for src,ix in base.groupby("source").groups.items():
        a=target.loc[ix];c=close.loc[ix]
        by_source.append({"source":src,"n_source":len(ix),"MIC_agreement_n":int(c.sum()),"joint_n":int(a.sum()),
                          "joint_fraction_of_source":a.mean(),"fraction_of_joint_group":a.sum()/target.sum(),
                          "joint_enrichment_vs_collection":a.mean()/target.mean()})
    save(by_source,"source_representation")
    check("source counts sum to subgroup size",sum(r["joint_n"] for r in by_source)==int(target.sum()))

    print("Comparing sequence profiles with matched controls",flush=True)
    # Two distinct controls: the full remaining collection, and the MUCH smaller
    # set that also agrees in MIC but does not cross the 50-point rank-gap rule.
    # Matching prediction bins conditions on model outputs; no causal claim.
    matches=[];support=[]
    for control_name,control in [("rest_of_collection",~target),("MIC_agreement_smaller_rank_gap",close&~target)]:
        check("target and control disjoint "+control_name,not (target&control).any())
        selection=target|control
        for bw in [None,1.,.5,.25]:
            dimensions=["source","length"]+(["median_MIC_bin"] if bw else [])
            bb=base.loc[selection].copy();bb["group"]=np.where(target.loc[selection],"target","control")
            if bw:bb["median_MIC_bin"]=np.floor(bb.median_log2_MIC/bw)
            count=bb.groupby(dimensions+["group"]).size().unstack("group").dropna()
            capacity=count.min(axis=1)
            support.append({"control":control_name,"matching":("source_length" if bw is None else "source_length_median_MIC"),
                            "log2_median_MIC_bin_width":bw,"target_n":int(target.sum()),"control_n":int(control.sum()),
                            "common_strata":len(count),"unweighted_matching_capacity":float(capacity.sum()),
                            "capacity_fraction_of_target":capacity.sum()/target.sum(),
                            "target_sequences_in_common_strata":int(count.target.sum()),
                            "fraction_target_in_common_strata":count.target.sum()/target.sum()})
            for weighting in ["each_sequence","equal_component_within_source"]:
                ww=pd.Series(1.,index=values.index) if weighting=="each_sequence" else 1/base.groupby(["source","component"]).source.transform("size")
                bb["weight"]=ww.loc[bb.index]
                sums=bb.groupby(dimensions+["group"]).weight.sum().unstack("group").dropna()
                mass=sums.min(axis=1)
                for feature in f:
                    bb["weighted_feature"]=f.loc[bb.index,feature].astype(float)*bb.weight
                    mean=(bb.groupby(dimensions+["group"]).weighted_feature.sum()/bb.groupby(dimensions+["group"]).weight.sum()).unstack("group")
                    left=mean.loc[mass.index,"target"];right=mean.loc[mass.index,"control"]
                    matches.append({"control":control_name,"log2_median_MIC_bin_width":bw,"weighting":weighting,"feature":feature,
                                    "common_strata":len(mass),"matched_mass":mass.sum(),"target_mean_on_common_strata":np.average(left,weights=mass),
                                    "control_mean_on_common_strata":np.average(right,weights=mass),
                                    "target_minus_control":np.average(left-right,weights=mass)})
    save(matches,"matched_sequence_profiles");save(support,"matched_comparison_support")
    mm=pd.DataFrame(matches)
    check("exact length matching gives zero length differences",np.allclose(mm[mm.feature.eq("length")].target_minus_control,0))

    # Population sensitivities retain the primary full-cohort rank reference, with
    # re-ranked populations separately labeled because the ranking question changes.
    corpus=set(pd.read_csv(ROOT/"scratchpad/replay/APEXDB.csv").Sequence)
    for line in (ROOT/"scratchpad/replay/ampredictor_train.txt").read_text().splitlines():
        fields=line.replace("\\t","\t").split("\t")
        if len(fields)>=2:corpus.add(fields[1].strip())
    order=np.argsort([hashlib.sha256(s.encode()).hexdigest() for s in values.index])
    ordered_components=component.iloc[order]
    representative=values.index[order[~ordered_components.duplicated().to_numpy()]]
    populations={"all":values.index,"cysteine_free_length8_to39":values.index[length.between(8,39)&f.contains_cysteine.eq(0)],
                 "exclude_exact_cached_training_corpora":values.index[~values.index.isin(corpus)],
                 "one_broader_component_representative":representative}
    sensitivities=[]
    for pop,idx in populations.items():
        for reference in ["original_full_cohort","recompute_within_population"]:
            gap=rank_gap.loc[idx] if reference=="original_full_cohort" else np.ptp(values.loc[idx].rank(pct=True),axis=1)
            near=close.loc[idx];joint=near&gap.ge(.5)
            sensitivities.append({"population":pop,"rank_reference":reference,"n":len(idx),"MIC_agreement_n":int(near.sum()),
                                  "joint_n":int(joint.sum()),"MIC_agreement_fraction":near.mean(),"joint_fraction":joint.mean(),
                                  "joint_fraction_of_MIC_agreement":joint.sum()/near.sum()})
    save(sensitivities,"population_sensitivity")
    cases=base.loc[close].join(f.drop(columns="length")).copy()
    cases["MIC_max_min_fold_ratio"]=np.exp2(spread.loc[close]);cases["rank_gap_percentage_points"]=100*rank_gap.loc[close]
    cases["large_rank_gap_subgroup"]=target.loc[close]
    for model in values:
        cases[model+"_MIC_uM"]=values.loc[close,model];cases[model+"_rank_percentile"]=100*ranks.loc[close,model]
    cases.to_csv(OUT/"peptides_with_similar_MIC.csv")
    selected=cases[cases.large_rank_gap_subgroup].sort_values("MIC_max_min_fold_ratio").drop_duplicates("component").head(12)
    selected.to_csv(OUT/"illustrative_examples.csv")
    target_components=component.loc[target].value_counts()
    residue_counts=Counter("".join(seq.loc[target]));residue_total=sum(residue_counts.values())
    save([{"AA":aa,"residue_count":residue_counts[aa],"fraction_of_all_subgroup_residues":residue_counts[aa]/residue_total} for aa in "ACDEFGHIKLMNPQRSTVWY"],"pooled_residue_composition")

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout="constrained")
    table=pd.DataFrame(prevalence).query("minimum_rank_gap_percentage_points == 50 and maximum_fold_MIC_spread >= 1.25 and maximum_fold_MIC_spread <= 8")
    xx=np.arange(len(table))
    axes[0].bar(xx-.18,table.MIC_agreement_fraction_of_all*100,.35,color="#b7c8d2",label="Similar MIC across all three")
    axes[0].bar(xx+.18,table.joint_fraction_of_all*100,.35,color="#225b78",label="Also rank gap ≥50 percentage points")
    axes[0].set(xticks=xx,xticklabels=table.maximum_fold_MIC_spread.map(lambda x:f"{x:g}×"),xlabel="Allowed largest/smallest MIC ratio",
                ylabel="Fraction of all 194,947 sequences (%)",title="A. Quantifying agreement and divergent rankings")
    axes[0].legend(fontsize=8)
    example=selected.iloc[0]
    positions=np.arange(3);models=list(values)
    axes[1].bar(positions,[example[m+"_rank_percentile"] for m in models],color=["#225b78","#c16936","#287a65"])
    for i,m in enumerate(models):axes[1].text(i,example[m+"_rank_percentile"]+2,f"{example[m+'_MIC_uM']:.2f} µM",ha="center",fontsize=9)
    axes[1].set(xticks=positions,xticklabels=models,ylim=(0,100),ylabel="Percentile in each model's list (lower = more favored)",
                title="B. Nearly identical MICs, very different priorities")
    fig.suptitle(str(selected.index[0])+": illustrative MIC agreement",fontsize=13,weight="bold")
    for ext in ["png","pdf"]:fig.savefig(OUT/("agreement_prevalence_and_example."+ext),dpi=185,bbox_inches="tight")
    plt.close(fig)
    save(CHECKS,"subanalysis_checks")
    summary={"input_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),"n_complete":n,
             "primary_max_MIC_fold_spread":2,"primary_min_rank_gap_percentage_points":50,
             "MIC_agreement_n":int(close.sum()),"MIC_agreement_fraction":float(close.mean()),"joint_n":int(target.sum()),
             "joint_fraction_of_all":float(target.mean()),"joint_fraction_of_MIC_agreement":float(target.sum()/close.sum()),
             "median_joint_rank_gap_percentage_points":float(100*rank_gap.loc[target].median()),
             "joint_group_component_n":len(target_components),"joint_group_top5_component_fraction":target_components.head(5).sum()/target.sum(),
             "median_joint_median_MIC_uM":float(np.exp2(z.loc[target].median(axis=1).median())),
             "median_other_agreement_median_MIC_uM":float(np.exp2(z.loc[close&~target].median(axis=1).median())),
             "all_joint_cases_APEX_has_lowest_rank_percentile":bool(ranks.loc[target].idxmin(axis=1).eq("APEX 1.1").all()),
             "subgroup_residue_total":residue_total,"L_K_R_fraction_of_pooled_residues":sum(residue_counts[a] for a in "LKR")/residue_total,
             "checks":len(CHECKS),"passed":sum(c["passed"] for c in CHECKS),"failed":[c for c in CHECKS if not c["passed"]],
             "python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"elapsed_seconds":time.monotonic()-start,
             "scope":"Exploratory joint subgroup defined from model outputs. No accuracy, unique motif, causal residue effect, or population significance inferred."}
    (OUT/"subanalysis_manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2));assert all(c["passed"] for c in CHECKS)


if __name__=="__main__":main()
