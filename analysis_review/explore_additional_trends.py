"""Exploratory trends beyond raw prediction-disagreement tails.

Run: python analysis_review/explore_additional_trends.py
All original inputs/outputs are preserved. New results: outputs/additional_trends.
No inference, downloads, significance claims, or experimental accuracy estimates.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import os
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/additional_trends"
OUT.mkdir(parents=True, exist_ok=True)
MAP = {"APEX 1.0":"apex_1_0_E._coli_ATCC11775", "APEX 1.1":"apex_1_1_E._coli_ATCC_11775",
       "MBC-Attention":"mbc_attention_E._coli_unspecified_strain", "AMPredictor":"ampredictor_E._coli_unspecified_strain"}
PRIMARY = ["APEX 1.1", "MBC-Attention", "AMPredictor"]
CLASS = ["macrel_general_AMP","ai4amp_general_AMP","amp_clip_general_AMP","amplify_general_AMP"]
CHECKS = []


def save(obj,name):
    pd.DataFrame(obj).to_csv(OUT/(name+".csv"),index=False)


def check(name,passed):
    CHECKS.append({"check":name,"passed":bool(passed)})


def weighted_rank(x,w):
    order=np.argsort(x,kind="stable");xx=x[order];ww=w[order]
    _,first,counts=np.unique(xx,return_index=True,return_counts=True)
    totals=np.add.reduceat(ww,first)
    middle=(np.cumsum(totals)-totals/2)/totals.sum()
    result=np.empty(len(x));result[order]=np.repeat(middle,counts)
    return result


def weighted_spearman(x,y,w):
    a=weighted_rank(np.asarray(x),np.asarray(w));b=weighted_rank(np.asarray(y),np.asarray(w))
    a=a-np.average(a,weights=w);b=b-np.average(b,weights=w)
    return np.average(a*b,weights=w)/np.sqrt(np.average(a*a,weights=w)*np.average(b*b,weights=w))


def figsave(fig,name):
    for ext in ["png","pdf"]:fig.savefig(OUT/(name+"."+ext),dpi=185,bbox_inches="tight")
    plt.close(fig)


def main():
    started=time.monotonic()
    path=ROOT/"df_damp_195625x57.tsv"
    raw=pd.read_csv(path,sep="\t").set_index("sequence",verify_integrity=True)
    source=pd.read_csv(ROOT/"scratchpad/source_sequences.tsv",sep="\t").set_index("sequence",verify_integrity=True).source
    component=pd.read_csv(ROOT/"analysis_review/outputs/sanity_checks/broader_sequence_components.csv").set_index("sequence",verify_integrity=True).component
    v=raw[list(MAP.values())].rename(columns={c:m for m,c in MAP.items()}).dropna().sort_index()
    z=np.log2(v);rank=z.rank(pct=True)
    f=pd.DataFrame({"source":source.reindex(z.index),"length":z.index.str.len(),"component":component.reindex(z.index)},index=z.index)
    f["signature"]=["".join(sorted(s)) for s in z.index]
    f["has_cysteine"]=z.index.str.contains("C")
    check("complete-case denominator",len(z)==194947)
    check("component coverage",f.component.notna().all())
    print("Exploring exact-composition sequence variants",flush=True)
    group=f.groupby(["source","signature"]).ngroup()
    counts=group.value_counts();same=counts[counts.ge(2)].index
    rows=[];details=[]
    for model in z:
        stats=z[model].groupby(group).agg(["min","max","size"])
        rr=rank[model].groupby(group).agg(["min","max"])
        for name,mask in [("all_same_source_composition_groups",stats["size"].ge(2)),
                          ("exactly_two_sequences_per_group",stats["size"].eq(2))]:
            b=stats.loc[mask];delta=b["max"]-b["min"]
            rows.append({"population":name,"model":model,"n_groups":len(b),"n_sequences":int(b["size"].sum()),
                         "median_log2_range":delta.median(),"median_fold_range":2**delta.median(),
                         "fraction_groups_ge_twofold_range":delta.ge(1).mean(),
                         "fraction_groups_ge_fourfold_range":delta.ge(2).mean(),
                         "median_global_percentile_range":(rr.loc[b.index,"max"]-rr.loc[b.index,"min"]).median()})
        for src,idx in f.groupby("source").groups.items():
            gids=np.intersect1d(group.loc[idx].unique(),same)
            b=stats.loc[gids];d=b["max"]-b["min"]
            details.append({"source":src,"model":model,"n_groups":len(b),"n_sequences":int(b["size"].sum()),
                            "median_log2_range":d.median(),"fraction_ge_twofold":d.ge(1).mean()})
    save(rows,"composition_matched_order_sensitivity");save(details,"composition_matched_by_source")
    # Fixed domain restriction, independently regrouping the restricted sequences.
    restricted=f.length.between(8,39)&~f.has_cysteine
    rgroup=f.loc[restricted].groupby(["source","signature"]).ngroup()
    domain=[]
    for model in z:
        b=z.loc[restricted,model].groupby(rgroup).agg(["min","max","size"]);b=b[b["size"].ge(2)]
        delta=b["max"]-b["min"]
        domain.append({"model":model,"n_groups":len(b),"n_sequences":int(b["size"].sum()),
                       "median_log2_range":delta.median(),"fraction_ge_twofold":delta.ge(1).mean()})
    save(domain,"composition_matched_domain_restriction")
    check("known exact-composition group count",int(counts.ge(2).sum())==6886)
    # Canonical orientation in two-member groups avoids choosing extreme examples
    # independently for every model. A model-independent pair set is shared by all.
    pair_rows=[];size_two=counts.index[counts.eq(2)]
    for gid,block in f.loc[group.isin(size_two)].groupby(group):
        left,right=sorted(block.index)
        check_counts=(f.at[left,"signature"]==f.at[right,"signature"] and f.at[left,"source"]==f.at[right,"source"])
        if not check_counts:raise AssertionError("Composition/source mismatch")
        row={"sequence_left":left,"sequence_right":right,"source":f.at[left,"source"],"length":len(left),
             "component_left":int(component[left]),"component_right":int(component[right])}
        for model in z:row[model+"_log2_right_minus_left"]=z.at[right,model]-z.at[left,model]
        pair_rows.append(row)
    pairs=pd.DataFrame(pair_rows);save(pairs,"exact_composition_two_member_pairs")
    check("two-member composition pair count",len(pairs)==5410)
    pair_agreement=[]
    for left,right in combinations(z,2):
        dl=pairs[left+"_log2_right_minus_left"];dr=pairs[right+"_log2_right_minus_left"]
        for margin in [0,.25,1.]:
            ok=dl.abs().gt(margin)&dr.abs().gt(margin)
            pair_agreement.append({"left":left,"right":right,"minimum_absolute_log2_change_each":margin,
                                   "n_pairs":len(pairs),"n_eligible":int(ok.sum()),
                                   "opposite_direction_fraction_among_eligible":dl.loc[ok].mul(dr.loc[ok]).lt(0).mean()})
    save(pair_agreement,"composition_matched_direction_agreement")

    # A position-specific subset: same observed core and added residue, present at
    # both termini. Deduplicate unordered longer-sequence pairs before reporting.
    seen_pairs=set();terminal=[];seqs=set(z.index)
    for core in z.index:
        if not 8<=len(core)<=49:continue
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            nn,cc=aa+core,core+aa
            if nn==cc or nn not in seqs or cc not in seqs or source[nn]!=source[cc]:continue
            key=tuple(sorted((nn,cc)))
            if key in seen_pairs:continue
            seen_pairs.add(key)
            rec={"core":core,"added_residue":aa,"N_extension":nn,"C_extension":cc,"source":source[nn]}
            for model in z:rec[model+"_log2_C_minus_N"]=z.at[cc,model]-z.at[nn,model]
            terminal.append(rec)
    term=pd.DataFrame(terminal);save(term,"terminal_position_composition_controls")
    terminal_summary=[]
    for model in z:
        delta=term[model+"_log2_C_minus_N"]
        terminal_summary.append({"model":model,"n_pairs":len(term),"median_signed_log2_C_minus_N":delta.median(),
                                 "median_absolute_log2_change":delta.abs().median(),"fraction_ge_twofold_change":delta.abs().ge(1).mean()})
    save(terminal_summary,"terminal_position_summary")

    print("Exploring classifier filtering of MIC-model favorites",flush=True)
    cp=raw.loc[z.index,CLASS]
    check("all classifiers available on MIC comparison cohort",cp.notna().all().all())
    filters=[];vote_rows=[]
    representative_order=np.argsort([hashlib.sha256(s.encode()).hexdigest() for s in z.index])
    ordered=f.iloc[representative_order];rep=ordered.index[~ordered.component.duplicated()]
    populations={"all":z.index,"cysteine_free_8_to_39":z.index[restricted],"one_broader_component_representative":rep}
    for pop,idx in populations.items():
        for threshold in [.3,.5,.7]:
            votes=cp.loc[idx].ge(threshold).sum(axis=1);allpos=votes.eq(4)
            k=int(np.ceil(len(idx)*.01));block=z.loc[idx]
            base=pd.DataFrame({"source":f.loc[idx,"source"],"length":f.loc[idx,"length"],"all_positive":allpos})
            baseline=base.groupby(["source","length"]).all_positive.transform("mean")
            for model in z:
                selected=block[model].nsmallest(k).index
                filters.append({"population":pop,"classifier_threshold":threshold,"model":model,"n":len(idx),"k":k,
                                "all_four_positive_cohort_fraction":allpos.mean(),"n_selected_retained_all_four":int(allpos.loc[selected].sum()),
                                "selected_retained_all_four_fraction":allpos.loc[selected].mean(),
                                "selected_all_four_negative_fraction":votes.loc[selected].eq(0).mean(),
                                "source_length_reference_positive_fraction":baseline.loc[selected].mean(),
                                "excess_above_source_length_reference":allpos.loc[selected].mean()-baseline.loc[selected].mean()})
            if pop=="all" and threshold==.5:
                pr=block[PRIMARY].rank(pct=True)
                for vv in range(5):
                    ix=votes.index[votes.eq(vv)]
                    vote_rows.append({"positive_classifier_votes":vv,"n":len(ix),
                                      "fraction_in_none_of_three_models_best10pct":pr.loc[ix].le(.1).sum(axis=1).eq(0).mean(),
                                      **{model+"_median_percentile":pr.loc[ix,model].median() for model in PRIMARY}})
    filter_table=pd.DataFrame(filters);save(filter_table,"classifier_filter_retention");save(vote_rows,"classifier_votes_vs_MIC_priority")

    print("Exploring matched-target activity-spectrum contrasts",flush=True)
    meta=pd.DataFrame(json.loads((ROOT/"scratchpad/damp_metadata.json").read_text())).set_index("column")
    older=meta[meta.model.eq("apex_1_0")]
    newer=meta[meta.model.eq("apex_1_1")&~meta.alias_unverified]
    old_cols=[];new_cols=[];names=[];strain=[]
    for col,row in newer.iterrows():
        old=older.index[older.strain_index.eq(row.strain_index)][0]
        if bool(older.loc[old,"alias_unverified"]):continue
        old_cols.append(old);new_cols.append(col);names.append(row.target);strain.append(int(row.strain_index))
    check("matched target panel has eight targets",len(names)==8)
    gn=np.isin(strain,[0,2,6,7]);gp=np.isin(strain,[3,26,27,28])
    check("matched Gram target counts",gn.sum()==gp.sum()==4 and np.all(gn|gp))
    z0=np.log2(raw[old_cols]);z1=np.log2(raw[new_cols]);z0.columns=names;z1.columns=names
    contrast={"APEX 1.0":z0.loc[:,gp].mean(axis=1)-z0.loc[:,gn].mean(axis=1),
              "APEX 1.1":z1.loc[:,gp].mean(axis=1)-z1.loc[:,gn].mean(axis=1)}
    target_pairs=[]
    for i,j in combinations(range(8),2):
        delta0=z0.iloc[:,j]-z0.iloc[:,i];delta1=z1.iloc[:,j]-z1.iloc[:,i]
        target_pairs.append({"target_left":names[i],"target_right":names[j],"n":len(raw),
                             "cross_version_spearman_of_log2_target_ratio":delta0.corr(delta1,method="spearman")})
    save(target_pairs,"target_pair_contrast_reproducibility")
    cross=[]
    for label,idx in [("all",raw.index)]+[(src,source.index[source.eq(src)]) for src in sorted(source.unique())]:
        g0=contrast["APEX 1.0"].loc[idx];g1=contrast["APEX 1.1"].loc[idx]
        eligible=g0.abs().gt(1)&g1.abs().gt(1)
        cross.append({"population":label,"n":len(idx),"cross_version_spearman":g0.corr(g1,method="spearman"),
                      "n_both_absolute_Gram_contrasts_gt_one_log2":int(eligible.sum()),
                      "same_direction_fraction_among_those":g0.loc[eligible].mul(g1.loc[eligible]).gt(0).mean()})
    save(cross,"Gram_contrast_reproducibility")
    omitted=[]
    for j,name in enumerate(names):
        gg=gn.copy();pp=gp.copy();gg[j]=False;pp[j]=False
        d0=z0.loc[:,pp].mean(axis=1)-z0.loc[:,gg].mean(axis=1)
        d1=z1.loc[:,pp].mean(axis=1)-z1.loc[:,gg].mean(axis=1)
        omitted.append({"omitted_target":name,"remaining_Gram_negative_targets":int(gg.sum()),
                        "remaining_Gram_positive_targets":int(pp.sum()),"cross_version_spearman":d0.corr(d1,method="spearman")})
    save(omitted,"Gram_contrast_leave_one_target_out")
    shifted0=z0+np.log2(512/140);shifted_contrast=shifted0.loc[:,gp].mean(axis=1)-shifted0.loc[:,gn].mean(axis=1)
    check("Gram contrast invariant to model-wide scale factor",np.allclose(shifted_contrast,contrast["APEX 1.0"]))
    row_shifted=z0.add(np.linspace(-3,3,len(z0)),axis=0)
    row_contrast=row_shifted.loc[:,gp].mean(axis=1)-row_shifted.loc[:,gn].mean(axis=1)
    check("Gram contrast removes arbitrary common per-peptide potency shifts",np.allclose(row_contrast,contrast["APEX 1.0"]))
    # A contrast cancels identical additive shifts across heads, but a rank-one
    # common axis with unequal target loadings can still generate that contrast.
    # Remove a cubic relationship between rank contrast and rank PC1 separately
    # in each version, then compare residual ranks by Pearson correlation.
    potency_rows=[];potency_residual=[]
    for version,zz in [("APEX 1.0",z0),("APEX 1.1",z1)]:
        pc=PCA(n_components=1,svd_solver="full").fit_transform(zz).ravel()
        pc_rank=pd.Series(pc,index=zz.index).rank(pct=True).to_numpy()
        contrast_rank=contrast[version].rank(pct=True).to_numpy()
        X=np.column_stack([np.ones(len(zz)),pc_rank,pc_rank**2,pc_rank**3])
        residual=contrast_rank-X@np.linalg.lstsq(X,contrast_rank,rcond=None)[0]
        potency_residual.append(residual)
        potency_rows.append({"version":version,"n":len(zz),
                             "Gram_contrast_vs_mean_log2_MIC_spearman":contrast[version].corr(zz.mean(axis=1),method="spearman"),
                             "absolute_Gram_contrast_vs_PC1_spearman":abs(np.corrcoef(contrast_rank,pc_rank)[0,1])})
    residual_correlation=float(np.corrcoef(potency_residual)[0,1])
    for row in potency_rows:row["cross_version_rank_residual_correlation_after_cubic_rank_PC1"]=residual_correlation
    save(potency_rows,"Gram_contrast_potency_axis_sensitivity")

    print("Checking charge associations under sequence/component aggregation",flush=True)
    seq=raw.index.to_series();length=seq.str.len()
    feat=pd.DataFrame({"source":source.reindex(raw.index),"component":component.reindex(raw.index),"length":length,
                       "charge_proxy":(seq.str.count("[KR]")-seq.str.count("[DE]"))/length,
                       "hydrophobic_fraction":seq.str.count("[AILMFWVY]")/length,"cysteine_fraction":seq.str.count("C")/length},index=raw.index)
    charge_rows=[];ven_sens=[];component_source_summary=[]
    for src,idx in feat.groupby("source").groups.items():
        b=feat.loc[idx].dropna(subset=["component"])
        counts_source=b.component.value_counts();w=1/b.component.map(counts_source).to_numpy()
        component_source_summary.append({"source":src,"n_sequences":len(b),"n_components":len(counts_source),
                                         "largest_component_n_in_source":int(counts_source.max()),
                                         "top10_components_fraction_of_source":counts_source.head(10).sum()/len(b)})
        for version,g in contrast.items():
            bb=b.assign(contrast=g.reindex(b.index))
            med=bb.groupby("component")[["charge_proxy","contrast"]].median()
            X=np.column_stack([np.ones(len(bb)),bb[["length","hydrophobic_fraction","cysteine_fraction"]].rank(pct=True)])
            Y=bb[["charge_proxy","contrast"]].rank(pct=True).to_numpy()
            residual=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0]
            ordinary=bb.charge_proxy.corr(bb.contrast,method="spearman")
            check("weighted Spearman recovers ordinary correlation "+src+"/"+version,
                  np.isclose(weighted_spearman(bb.charge_proxy,bb.contrast,np.ones(len(bb))),ordinary))
            charge_rows.append({"source":src,"version":version,"n_sequences":len(bb),"n_components":len(med),
                                "row_spearman_charge_vs_Gram_contrast":ordinary,
                                "rank_residual_correlation_after_length_hydrophobicity_cysteine":np.corrcoef(residual.T)[0,1],
                                "equal_component_weighted_spearman":weighted_spearman(bb.charge_proxy,bb.contrast,w),
                                "component_median_spearman":med.charge_proxy.corr(med.contrast,method="spearman")})
            if src=="Venomics":
                for drop in [0,1,3,5,10]:
                    keep=~bb.component.isin(counts_source.head(drop).index) if drop else pd.Series(True,index=bb.index)
                    ven_sens.append({"version":version,"largest_components_removed":drop,"remaining_sequences":int(keep.sum()),
                                     "charge_vs_Gram_contrast_spearman":bb.loc[keep,"charge_proxy"].corr(bb.loc[keep,"contrast"],method="spearman")})
    save(charge_rows,"charge_Gram_contrast_aggregation_sensitivity");save(ven_sens,"Venomics_large_component_sensitivity")
    save(component_source_summary,"source_sequence_component_concentration")
    # Compare Gram contrasts at component medians on the same complete-case set.
    cb=pd.DataFrame({version:g.reindex(z.index) for version,g in contrast.items()}).groupby(component.reindex(z.index)).median()
    save([{"population":"complete_case_component_medians","n_components":len(cb),
           "cross_version_Gram_contrast_spearman":cb.iloc[:,0].corr(cb.iloc[:,1],method="spearman")}],"Gram_contrast_component_reproducibility")

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    colors=["#9c7770","#225b78","#c16936","#287a65"]
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    order=pd.DataFrame(rows);b=order[order.population.eq("exactly_two_sequences_per_group")].set_index("model").loc[list(z)]
    axes[0].barh(b.index,b.fraction_groups_ge_twofold_range*100,color=colors)
    axes[0].set(xlabel="Composition-matched pairs with ≥2-fold prediction change (%)",title="A. Sequence order matters differently across models");axes[0].invert_yaxis()
    for i,r in enumerate(b.itertuples()):axes[0].text(r.fraction_groups_ge_twofold_range*100+.3,i,f"{r.fraction_groups_ge_twofold_range:.1%}",va="center")
    axes[0].set_xlim(0,25)
    b=filter_table[filter_table.population.eq("all")&filter_table.classifier_threshold.eq(.5)].set_index("model").loc[list(z)]
    axes[1].barh(b.index,b.selected_retained_all_four_fraction*100,color=colors)
    axes[1].scatter(b.source_length_reference_positive_fraction*100,np.arange(4),marker="|",s=140,color="black",label="Source/length reference")
    axes[1].set(xlim=(0,100),xlabel="Original lowest-MIC 1% list passing all four classifiers (%)",title="B. A classifier-consensus filter favors APEX-selected lists");axes[1].invert_yaxis();axes[1].legend(fontsize=8)
    fig.suptitle("Composition alone and classifier consensus do not define a common MIC ranking",fontsize=13,weight="bold")
    figsave(fig,"figure_7_sequence_order_and_filtering")
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    hx=axes[0].hexbin(contrast["APEX 1.0"],contrast["APEX 1.1"],gridsize=65,bins="log",cmap="Blues",mincnt=1)
    axes[0].axhline(0,color="gray",lw=.7);axes[0].axvline(0,color="gray",lw=.7)
    axes[0].set(xlabel="APEX 1.0 mean log2 MIC(GP) − mean log2 MIC(GN)",ylabel="APEX 1.1: same target contrast",
                title=f"A. Relative target preferences persist (ρ={cross[0]['cross_version_spearman']:.3f})")
    fig.colorbar(hx,ax=axes[0],label="Sequences per bin",shrink=.75)
    q=pd.DataFrame(charge_rows).query("source == 'Venomics'")
    labels=["Each sequence equal","Each component equal","Medians per component"]
    for i,(version,color) in enumerate([("APEX 1.0","#9c7770"),("APEX 1.1","#225b78")]):
        r=q[q.version.eq(version)].iloc[0]
        y=[r.row_spearman_charge_vs_Gram_contrast,r.equal_component_weighted_spearman,r.component_median_spearman]
        axes[1].bar(np.arange(3)+(i-.5)*.3,y,width=.29,color=color,label=version)
    axes[1].axhline(0,color="gray",lw=.8)
    axes[1].set(xticks=np.arange(3),xticklabels=labels,ylabel="Charge-proxy / Gram-contrast rank association",ylim=(-.35,.4),
                title="B. Venomics charge association changes with aggregation")
    axes[1].tick_params(axis="x",labelrotation=15);axes[1].legend(fontsize=8)
    fig.suptitle("Target contrasts reveal shared model structure and dependence-sensitive associations",fontsize=13,weight="bold")
    figsave(fig,"figure_8_spectrum_and_Venomics_aggregation")
    save(CHECKS,"additional_trend_checks")
    summary={"input_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"complete_cases":len(z),
             "same_source_composition_groups":int(counts.ge(2).sum()),"sequences_in_those_groups":int(counts[counts.ge(2)].sum()),
             "two_member_groups":len(pairs),"unique_terminal_position_pairs":len(term),
             "Gram_contrast_cross_version_rho":cross[0]["cross_version_spearman"],
             "Gram_contrast_cross_version_rank_residual_correlation_after_PC1":residual_correlation,
             "Gram_contrast_leave_target_out_rho_min":min(r["cross_version_spearman"] for r in omitted),
             "Gram_contrast_leave_target_out_rho_max":max(r["cross_version_spearman"] for r in omitted),
             "checks":len(CHECKS),"passed":sum(c["passed"] for c in CHECKS),"failed":[c for c in CHECKS if not c["passed"]],
             "elapsed_seconds":time.monotonic()-started,
             "dependency_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
                 ROOT/"scratchpad/source_sequences.tsv",ROOT/"scratchpad/damp_metadata.json",
                 ROOT/"analysis_review/outputs/sanity_checks/broader_sequence_components.csv"]},
             "scope":"Exploratory selected-collection model behavior. No accuracy, causal residue effects, or population significance inferred."}
    (OUT/"additional_trends_manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2));print(pd.DataFrame(rows).to_string(index=False));print(q.to_string(index=False))
    assert all(c["passed"] for c in CHECKS)


if __name__=="__main__":main()
