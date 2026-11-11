"""Sequence-local analysis of disagreement in the existing MIC predictions.

Run: python analysis_review/sequence_disagreement.py
No model inference or network. New results go to outputs/sequence_disagreement.
All original files and previous outputs are preserved.
"""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
from collections import Counter
import hashlib
import json
import os
import platform
import re
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"venomics-mpl"))
os.environ["MPLBACKEND"]="Agg"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"analysis_review/outputs/sequence_disagreement"
OUT.mkdir(parents=True,exist_ok=True)
MAP={"APEX 1.1":"apex_1_1_E._coli_ATCC_11775","MBC-Attention":"mbc_attention_E._coli_unspecified_strain",
     "AMPredictor":"ampredictor_E._coli_unspecified_strain"}
MODELS=list(MAP)
CHECKS=[]


def save(obj,name):
    pd.DataFrame(obj).to_csv(OUT/(name+".csv"),index=False)


def check(name,condition):
    CHECKS.append({"check":name,"passed":bool(condition)})


def fold_of_component(c):
    return hashlib.sha256(("sequence-disagreement-20261106:"+str(int(c))).encode()).digest()[0]%2


def absolute_distance_sums(x):
    ordered=np.sort(x);prefix=np.r_[0,np.cumsum(ordered)]
    j=np.searchsorted(ordered,x,side="left")
    return x*j-prefix[j]+prefix[-1]-prefix[j]-x*(len(x)-j)


def grouped_distance_sums(x,codes):
    """Exact within-group sums of absolute scalar distances without pair matrices."""
    x=np.asarray(x);codes=np.asarray(codes)
    order=np.lexsort((x,codes));xx=x[order];gg=codes[order]
    first=np.r_[0,np.flatnonzero(gg[1:]!=gg[:-1])+1]
    count=np.diff(np.r_[first,len(xx)])
    prefix=np.r_[0,np.cumsum(xx)]
    lower_sum=prefix[:-1]-np.repeat(prefix[first],count)
    within_position=np.arange(len(xx))-np.repeat(first,count)
    totals=np.add.reduceat(xx,first)
    sorted_sum=xx*(2*within_position-np.repeat(count,count))+np.repeat(totals,count)-2*lower_sum
    result=np.empty(len(x));result[order]=np.maximum(sorted_sum,0.)
    return result


@lru_cache(maxsize=None)
def order_features(s):
    def longest(pattern,string):
        return max([len(x) for x in re.findall(pattern,string)]+[0])
    return np.array([
        longest("[KR]+",s),longest("[AILMFWVY]+",s),
        sum(s[j] in "KR" and s[j+1] in "KR" for j in range(len(s)-1)),
        sum(c in "KR" for c in s[:3])-sum(c in "KR" for c in s[-3:]),
        min(len(s),longest("[AILMFWVY]+",s+s)),
    ],dtype=float)


def main():
    started=time.monotonic();path=ROOT/"df_damp_195625x57.tsv"
    raw=pd.read_csv(path,sep="\t",usecols=["sequence"]+list(MAP.values())).set_index("sequence",verify_integrity=True)
    values=raw.rename(columns={v:k for k,v in MAP.items()}).dropna().sort_index()
    z=np.log2(values);ranks=z.rank(pct=True);centered=z-z.median()
    src=pd.read_csv(ROOT/"scratchpad/source_sequences.tsv",sep="\t").set_index("sequence",verify_integrity=True).source.reindex(z.index)
    comp=pd.read_csv(ROOT/"analysis_review/outputs/sanity_checks/broader_sequence_components.csv").set_index("sequence",verify_integrity=True).component.reindex(z.index)
    f=pd.DataFrame({"source":src,"length":z.index.str.len(),"component":comp},index=z.index)
    scores=pd.DataFrame({"raw_log2_range":np.ptp(z,axis=1),"centered_log2_range":np.ptp(centered,axis=1),
                         "rank_range":np.ptp(ranks,axis=1)},index=z.index)
    check("primary comparison has 194947 finite positive sequences",len(z)==194947 and np.isfinite(z).all().all())
    check("all sequences mapped to source and observed component",f.notna().all().all())
    print("Enumerating exact one-substitution pairs within sources",flush=True)
    by_length={}
    for s in z.index:
        if len(s)>=8:by_length.setdefault(len(s),[]).append(s)
    records=[]
    for length,sequences in by_length.items():
        for pos in range(length):
            buckets={}
            for s in sequences:buckets.setdefault(s[:pos]+s[pos+1:],[]).append(s)
            for members in buckets.values():
                if len(members)<2:continue
                for left,right in combinations(members,2):
                    if src[left]!=src[right]:continue
                    if left[pos]>right[pos]:left,right=right,left
                    records.append((left,right,pos+1,length,left[pos],right[pos],src[left],int(comp[left])))
    pairs=pd.DataFrame(records,columns=["sequence_from","sequence_to","position_1based","length","AA_from","AA_to","source","component"])
    check("all substitution partners share their observed component",bool(np.array_equal(comp.loc[pairs.sequence_to].to_numpy(),pairs.component.to_numpy())))
    check("no duplicate substitution pairs",not pairs.duplicated(["sequence_from","sequence_to"]).any())
    check("exact mutation count from preliminary enumeration",len(pairs)==41253)
    for row in pairs[["sequence_from","sequence_to","position_1based"]].itertuples(index=False):
        s,t=row.sequence_from,row.sequence_to
        if sum(a!=b for a,b in zip(s,t))!=1 or len(s)!=len(t):raise AssertionError("not one substitution")
        if s[row.position_1based-1]==t[row.position_1based-1]:raise AssertionError("wrong mutation site")
    check("all one-substitution identities and sites verified directly",True)
    pairs["mutation"]=pairs.AA_from+"→"+pairs.AA_to
    pairs["component_subset"]=pairs.component.map(fold_of_component)
    for model in MODELS:
        pairs[model+"_delta_log2"]=z.loc[pairs.sequence_to,model].to_numpy()-z.loc[pairs.sequence_from,model].to_numpy()
    for score in scores:
        pairs[score+"_from"]=scores.loc[pairs.sequence_from,score].to_numpy()
        pairs[score+"_to"]=scores.loc[pairs.sequence_to,score].to_numpy()
        pairs[score+"_delta"]=pairs[score+"_to"]-pairs[score+"_from"]
    for left,right in combinations(MODELS,2):
        name=left+" minus "+right+"_gap_delta_log2"
        pairs[name]=pairs[left+"_delta_log2"]-pairs[right+"_delta_log2"]
        baseline=z[left]-z[right]
        independent=baseline.loc[pairs.sequence_to].to_numpy()-baseline.loc[pairs.sequence_from].to_numpy()
        check("mutation gap identity "+left+"/"+right,np.allclose(pairs[name],independent))
        centered_gap=centered[left]-centered[right]
        check("gap changes unaffected by model-wise offsets "+left+"/"+right,
              np.allclose(independent,centered_gap.loc[pairs.sequence_to].to_numpy()-centered_gap.loc[pairs.sequence_from].to_numpy()))

    # A fixed-anchor reference: expected distance to a uniformly selected OTHER
    # sequence of the same source and length. The expectation is exact, not a
    # Monte Carlo null. It includes related alternatives, making no independence
    # claim, and a second reference excludes the anchor's entire component.
    print("Checking local clustering against source/length references",flush=True)
    references=pd.DataFrame(index=z.index)
    outer_codes=f.groupby(["source","length"]).ngroup().to_numpy()
    inner_codes=f.groupby(["source","length","component"]).ngroup().to_numpy()
    n_total=f.groupby(["source","length"]).source.transform("size")
    n_within=f.groupby(["source","length","component"]).source.transform("size")
    for score in scores:
        total_sums=pd.Series(grouped_distance_sums(scores[score].to_numpy(),outer_codes),index=z.index)
        within_sums=pd.Series(grouped_distance_sums(scores[score].to_numpy(),inner_codes),index=z.index)
        references[score+"_all"]=total_sums/(n_total-1).replace(0,np.nan)
        references[score+"_outside_component"]=(total_sums-within_sums)/(n_total-n_within).replace(0,np.nan)
    for score in scores:
        pairs[score+"_absolute_change"]=pairs[score+"_delta"].abs()
        for ref in ["all","outside_component"]:
            pairs[score+"_expected_"+ref]=(references.loc[pairs.sequence_from,score+"_"+ref].to_numpy()+references.loc[pairs.sequence_to,score+"_"+ref].to_numpy())/2
    tiny=scores.iloc[:10,2].to_numpy();dist=abs(tiny[:,None]-tiny[None,:]).sum(axis=1)
    check("absolute-distance formula equals explicit matrix",np.allclose(absolute_distance_sums(tiny),dist))
    toy=np.array([.3,.1,.1,.7,.2,.8,.3]);toy_codes=np.array([2,1,2,1,2,1,3])
    brute=np.array([sum(abs(x-toy[toy_codes==c])) for x,c in zip(toy,toy_codes)])
    check("grouped distance formula equals explicit groups with ties",np.allclose(grouped_distance_sums(toy,toy_codes),brute))
    locality=[]
    for score in scores:
        for name,b in [("all",pairs)]+list(pairs.groupby("source")):
            for ref in ["all","outside_component"]:
                valid=b.dropna(subset=[score+"_expected_"+ref])
                for weighting in ["each_pair","equal_component"]:
                    actual=valid[score+"_absolute_change"];expected=valid[score+"_expected_"+ref]
                    if weighting=="equal_component":
                        tmp=valid.groupby("component")[[score+"_absolute_change",score+"_expected_"+ref]].mean()
                        observed_mean=tmp.iloc[:,0].mean();reference_mean=tmp.iloc[:,1].mean()
                    else:observed_mean=actual.mean();reference_mean=expected.mean()
                    locality.append({"population":name,"score":score,"reference":ref,"weighting":weighting,
                                     "n_pairs":len(valid),"n_components":valid.component.nunique(),
                                     "mean_observed_absolute_change":observed_mean,"mean_reference_absolute_change":reference_mean,
                                     "observed_reference_ratio":observed_mean/reference_mean})
    save(locality,"local_disagreement_coherence")
    opposite=[]
    for left,right in combinations(MODELS,2):
        dl=pairs[left+"_delta_log2"];dr=pairs[right+"_delta_log2"]
        for margin in [.5,1.]:
            eligible=dl.abs().ge(margin)&dr.abs().ge(margin);opposed=eligible&dl.mul(dr).lt(0)
            ff=pd.DataFrame({"component":pairs.component,"opposed":opposed,"eligible":eligible})
            opposite.append({"left":left,"right":right,"minimum_absolute_log2_response_each":margin,
                             "n_pairs":len(pairs),"n_components":pairs.component.nunique(),"n_both_responses_above_margin":int(eligible.sum()),
                             "n_strong_opposite":int(opposed.sum()),"fraction_all_pairs_strong_opposite":opposed.mean(),
                             "fraction_eligible_pairs_opposite":opposed.loc[eligible].mean(),
                             "n_components_with_any_opposite":pairs.loc[opposed,"component"].nunique(),
                             "equal_component_mean_fraction_all_pairs_opposite":ff.groupby("component").opposed.mean().mean()})
    save(opposite,"single_substitution_opposite_responses")
    global_change=[]
    for score,threshold in [("raw_log2_range",1.),("centered_log2_range",1.),("rank_range",.1)]:
        indicator=pairs[score+"_delta"].abs().ge(threshold)
        global_change.append({"score":score,"change_threshold":threshold,"n_pairs":len(pairs),"n_at_least_threshold":int(indicator.sum()),
                              "fraction_pairs":indicator.mean(),"equal_component_fraction":pd.DataFrame({"c":pairs.component,"v":indicator}).groupby("c").v.mean().mean()})
    save(global_change,"single_substitution_score_changes")
    # Report the complete mutation-type screen; no p-values or unreported search
    # for a preferred residue. Small cells remain explicitly visible.
    mutation_summary=[];mutation_source=[]
    for mutation,b in pairs.groupby("mutation"):
        row={"mutation":mutation,"n_pairs":len(b),"n_components":b.component.nunique(),"n_sources":b.source.nunique()}
        for model in MODELS:
            row[model+"_median_response_log2"]=b[model+"_delta_log2"].median()
            row[model+"_equal_component_mean_response_log2"]=b.groupby("component")[model+"_delta_log2"].mean().mean()
        for score in scores:
            row[score+"_median_delta"]=b[score+"_delta"].median()
            row[score+"_equal_component_mean_delta"]=b.groupby("component")[score+"_delta"].mean().mean()
        mutation_summary.append(row)
        for (src,subset),bb in b.groupby(["source","component_subset"]):
            mutation_source.append({"mutation":mutation,"source":src,"component_subset":int(subset),"n_pairs":len(bb),
                                    "n_components":bb.component.nunique(),"rank_range_mean_delta":bb.rank_range_delta.mean(),
                                    "centered_log2_range_mean_delta":bb.centered_log2_range_delta.mean()})
    save(mutation_summary,"all_mutation_type_summaries");save(mutation_source,"mutation_source_and_component_subsets")
    # First-cysteine introduction: an observed edit, not a test of model validity
    # outside a documented training composition.
    cys=[]
    for row in pairs.itertuples(index=False):
        if row.AA_from=="C" or row.AA_to=="C":
            first=row.sequence_to if row.AA_from=="C" else row.sequence_from
            if "C" not in first:cys.append(row.Index if hasattr(row,"Index") else (row.sequence_from,row.sequence_to))
    first_cys=pairs.set_index(["sequence_from","sequence_to"]).loc[cys].reset_index() if cys else pairs.iloc[:0].copy()
    cys_summary=[];cys_source=[]
    if len(first_cys):
        direction=np.where(first_cys.AA_to.eq("C"),1.,-1.)
        for score in scores:
            delta=first_cys[score+"_delta"].to_numpy()*direction
            cb=first_cys.assign(oriented_delta=delta)
            cys_summary.append({"score":score,"n_pairs":len(first_cys),"n_components":first_cys.component.nunique(),
                                "median_first_C_introduction_delta":float(np.median(delta)),"mean_delta":float(np.mean(delta)),
                                "equal_component_mean_delta":cb.groupby("component").oriented_delta.mean().mean(),
                                "top_three_components_fraction_of_pairs":cb.component.value_counts().head(3).sum()/len(cb)})
            for (src,subset),b in cb.groupby(["source","component_subset"]):
                cys_source.append({"score":score,"source":src,"component_subset":int(subset),"n_pairs":len(b),
                                   "n_components":b.component.nunique(),"mean_delta":b.oriented_delta.mean(),
                                   "equal_component_mean_delta":b.groupby("component").oriented_delta.mean().mean()})
    save(cys_summary,"first_cysteine_introduction")
    save(cys_source,"first_cysteine_source_and_component_subsets")
    pairs.to_csv(OUT/"single_substitution_pairs.csv.gz",index=False,compression={"method":"gzip","mtime":0})

    print("Testing sequence arrangement with composition held fixed",flush=True)
    order=pd.read_csv(ROOT/"analysis_review/outputs/additional_trends/exact_composition_two_member_pairs.csv")
    parent={int(c):int(c) for c in set(order.component_left)|set(order.component_right)}
    def root(c):
        c=int(c)
        while parent[c]!=c:parent[c]=parent[parent[c]];c=parent[c]
        return c
    for a,b in zip(order.component_left,order.component_right):parent[root(a)]=root(b)
    order["pair_component"]=order.component_left.map(root)
    order["component_subset"]=order.pair_component.map(fold_of_component)
    kinds=[];all_adjacent=[]
    for s,t in zip(order.sequence_left,order.sequence_right):
        changed=[i for i,(a,b) in enumerate(zip(s,t)) if a!=b]
        all_adjacent.append(len(changed)==2 and changed[1]-changed[0]==1)
        if t in s+s:kind="letter_rotation"
        elif len(changed)==2 and changed[1]-changed[0]==1:kind="adjacent_swap"
        elif len(changed)==2:kind="separated_swap"
        else:kind="other_reordering"
        kinds.append(kind)
    order["edit_kind"]=kinds
    for score in scores:
        order[score+"_delta"]=scores.loc[order.sequence_right,score].to_numpy()-scores.loc[order.sequence_left,score].to_numpy()
    af=np.array([order_features(s) for s in order.sequence_left]);bf=np.array([order_features(s) for s in order.sequence_right]);difference=bf-af
    rotations=order.edit_kind.eq("letter_rotation").to_numpy()
    check("circular run control invariant to letter rotations",np.all(difference[rotations,4]==0))
    edit_summary=[]
    for kind,b in order.groupby("edit_kind"):
        for model in MODELS:
            delta=b[model+"_log2_right_minus_left"]
            edit_summary.append({"edit_kind":kind,"model":model,"n_pairs":len(b),"n_pair_components":b.pair_component.nunique(),
                                 "median_absolute_log2_response":delta.abs().median(),"n_response_at_least_twofold":int(delta.abs().ge(1).sum()),
                                 "fraction_response_at_least_twofold":delta.abs().ge(1).mean()})
    save(edit_summary,"composition_matched_edit_classes")
    adjacent_summary=[]
    for model in MODELS:
        delta=order.loc[all_adjacent,model+"_log2_right_minus_left"]
        adjacent_summary.append({"model":model,"n_all_adjacent_swap_pairs":len(delta),
                                 "n_also_classified_as_letter_rotation":int(np.sum(np.asarray(all_adjacent)&rotations)),
                                 "n_response_at_least_twofold":int(delta.abs().ge(1).sum()),
                                 "fraction_response_at_least_twofold":delta.abs().ge(1).mean()})
    save(adjacent_summary,"all_adjacent_swap_responses")
    pattern_summary=[];pattern_source=[];pattern_subsets=[]
    feature_names=["longest_basic_run","longest_hydrophobic_run","adjacent_basic_pairs","basic_residues_first3_minus_last3","circular_hydrophobic_run_control"]
    order_domain=(order.length.between(8,39)&~order.sequence_left.str.contains("C")).to_numpy()
    for j,feature in enumerate(feature_names):
        for category,mask in [("all",np.ones(len(order),dtype=bool)),("not_letter_rotations",~rotations),
                              ("cysteine_free_8_to_39",order_domain),("cysteine_free_non_rotations",order_domain&~rotations)]:
            valid=mask&(difference[:,j]!=0);bb=order.loc[valid].copy();direction=np.sign(difference[valid,j])
            for score in scores:
                bb["oriented_delta"]=bb[score+"_delta"].to_numpy()*direction
                pattern_summary.append({"feature":feature,"pair_selection":category,"score":score,"n_pairs":len(bb),
                                        "n_pair_components":bb.pair_component.nunique(),"median_more_minus_less":bb.oriented_delta.median(),
                                        "mean_more_minus_less":bb.oriented_delta.mean(),
                                        "equal_pair_component_mean_more_minus_less":bb.groupby("pair_component").oriented_delta.mean().mean()})
                for subset,b in bb.groupby("component_subset"):
                    pattern_subsets.append({"feature":feature,"pair_selection":category,"score":score,"component_subset":int(subset),
                                            "n_pairs":len(b),"n_pair_components":b.pair_component.nunique(),
                                            "median_more_minus_less":b.oriented_delta.median(),"mean_more_minus_less":b.oriented_delta.mean(),
                                            "equal_component_mean_more_minus_less":b.groupby("pair_component").oriented_delta.mean().mean()})
                if score=="rank_range":
                    for (src,subset),b in bb.groupby(["source","component_subset"]):
                        pattern_source.append({"feature":feature,"pair_selection":category,"source":src,"component_subset":int(subset),
                                               "n_pairs":len(b),"n_pair_components":b.pair_component.nunique(),
                                               "median_more_minus_less_rank_range":b.oriented_delta.median(),
                                               "mean_more_minus_less_rank_range":b.oriented_delta.mean(),
                                               "equal_component_mean_rank_range_delta":b.groupby("pair_component").oriented_delta.mean().mean()})
    save(pattern_summary,"composition_matched_arrangement_features");save(pattern_source,"arrangement_source_and_component_subsets")
    save(pattern_subsets,"arrangement_disjoint_component_subsets")
    order[["sequence_left","sequence_right","source","length","pair_component","edit_kind"]].to_csv(OUT/"order_pair_classification.csv",index=False)

    print("Comparing disagreement profiles within source and length",flush=True)
    seq=z.index.to_series();length=f.length
    features=pd.DataFrame({"basic_fraction":seq.str.count("[KR]")/length,"acidic_fraction":seq.str.count("[DE]")/length,
                           "charge_proxy":(seq.str.count("[KR]")-seq.str.count("[DE]"))/length,
                           "hydrophobic_fraction":seq.str.count("[AILMFWVY]")/length,"cysteine_fraction":seq.str.count("C")/length,
                           "proline_fraction":seq.str.count("P")/length,"glycine_fraction":seq.str.count("G")/length},index=z.index)
    features["composition_entropy_bits"]=[-sum((n/len(s))*np.log2(n/len(s)) for n in Counter(s).values()) for s in z.index]
    profiles=[]
    for score in scores:
        low,high=scores[score].quantile([.05,.95])
        groups=pd.Series(np.select([scores[score].le(low),scores[score].ge(high)],["low","high"],default="middle"),index=z.index)
        for weighting in ["each_sequence","equal_component_within_source"]:
            weights=pd.Series(1.,index=z.index) if weighting=="each_sequence" else 1/f.groupby(["source","component"]).source.transform("size")
            b=f.join(features).assign(group=groups,weight=weights);b=b[b.group.ne("middle")].copy()
            for feature in features:
                b["weighted_feature"]=b[feature]*b.weight
                aggregate=b.groupby(["source","length","group"])[["weight","weighted_feature"]].sum()
                mean=(aggregate.weighted_feature/aggregate.weight).unstack("group")
                weight=aggregate.weight.unstack("group").dropna();capacity=weight.min(axis=1)
                delta=mean.loc[capacity.index,"high"]-mean.loc[capacity.index,"low"]
                for name,mask in [("all",np.ones(len(capacity),dtype=bool))]+[(src,capacity.index.get_level_values("source")==src) for src in sorted(f.source.unique())]:
                    if not np.any(mask):continue
                    profiles.append({"score":score,"weighting":weighting,"population":name,"feature":feature,
                                     "n_common_source_length_strata":int(np.sum(mask)),"matched_weight_capacity":float(capacity.loc[mask].sum()),
                                     "mean_high_minus_low_in_common_strata":float(np.average(delta.loc[mask],weights=capacity.loc[mask]))})
    save(profiles,"source_length_matched_sequence_profiles")
    # Concrete examples, one per observed component; include a comprehensible
    # fixed example, then distinct-component extremes in a conservative domain.
    eligible=pairs[pairs.sequence_from.str.count("C").eq(0)&pairs.sequence_to.str.count("C").eq(0)&pairs.length.between(8,39)].copy()
    eligible["selection_magnitude"]=(eligible["APEX 1.1_delta_log2"]-eligible["AMPredictor_delta_log2"]).abs()
    mandatory=pd.concat([
        eligible[(eligible.sequence_from=="SLLDLIRKLLT")&(eligible.sequence_to=="SLLILIRKLLT")],
        eligible[(eligible.sequence_from=="LLKKLLKKLLKK")&(eligible.sequence_to=="LLKKLLKKLLKL")],
    ])
    candidate=pd.concat([mandatory,eligible.sort_values("selection_magnitude",ascending=False)]).drop_duplicates("component").head(8)
    training={}
    training["APEXDB"]=set(pd.read_csv(ROOT/"scratchpad/replay/APEXDB.csv").Sequence)
    amp=set()
    for line in (ROOT/"scratchpad/replay/ampredictor_train.txt").read_text().splitlines():
        fields=line.replace("\\t","\t").split("\t")
        if len(fields)>=2:amp.add(fields[1].strip())
    training["AMPredictor_train"]=amp
    examples=[]
    for _,row in candidate.iterrows():
        rec={"sequence_from":row.sequence_from,"sequence_to":row.sequence_to,"edit":f"{row.AA_from}{row.position_1based}{row.AA_to}",
             "source":row.source,"component":int(row.component),"length":int(row.length)}
        for model in MODELS:
            rec[model+"_from_uM"]=values.at[row.sequence_from,model];rec[model+"_to_uM"]=values.at[row.sequence_to,model]
        for corpus,ss in training.items():
            rec[corpus+"_from_exact_overlap"]=row.sequence_from in ss;rec[corpus+"_to_exact_overlap"]=row.sequence_to in ss
        examples.append(rec)
    save(examples,"single_substitution_examples")

    print("Writing sequence-level figures",flush=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    loc=pd.DataFrame(locality);loc=loc[(loc.score=="rank_range")&(loc.reference=="all")&(loc.weighting=="each_pair")&(loc.population!="all")]
    axes[0].barh(loc.population.str.replace("_"," "),loc.mean_reference_absolute_change*100,color="#d8e1e6",label="Same source and length reference")
    axes[0].barh(loc.population.str.replace("_"," "),loc.mean_observed_absolute_change*100,color="#225b78",label="One-letter neighbors")
    axes[0].set(xlabel="Mean difference in rank-disagreement score (percentage points)",title="A. Disagreement is locally similar in sequence space");axes[0].legend(fontsize=7.5)
    ex=pd.DataFrame(examples).iloc[0]
    for j,(model,color) in enumerate(zip(MODELS,["#225b78","#c16936","#287a65"])):
        axes[1].plot([0,1],[ex[model+"_from_uM"],ex[model+"_to_uM"]],"o-",label=model,color=color,lw=2)
    axes[1].set(yscale="log",xticks=[0,1],xticklabels=[ex.sequence_from,ex.sequence_to],ylabel="Predicted MIC (µM; lower means more potent)",
                title=f"B. A single {ex.edit[0]}→{ex.edit[-1]} edit can drive opposing predictions")
    axes[1].tick_params(axis="x",labelrotation=15);axes[1].legend(fontsize=8)
    for ext in ["png","pdf"]:fig.savefig(OUT/("sequence_locality_and_edit_example."+ext),dpi=185,bbox_inches="tight")
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout="constrained")
    es=pd.DataFrame(edit_summary);pivot=es.pivot(index="edit_kind",columns="model",values="fraction_response_at_least_twofold").reindex(["adjacent_swap","separated_swap","letter_rotation","other_reordering"])
    pivot.mul(100).plot.barh(ax=axes[0],color={"APEX 1.1":"#225b78","MBC-Attention":"#c16936","AMPredictor":"#287a65"})
    axes[0].set(xlabel="Pairs with ≥2-fold predicted MIC response (%)",ylabel="",title="A. Equal composition: edit type matters")
    arr=pd.DataFrame(pattern_summary);arr=arr[(arr.score=="rank_range")&arr.pair_selection.eq("all")&~arr.feature.eq("circular_hydrophobic_run_control")]
    axes[1].barh(arr.feature.str.replace("_"," "),arr.mean_more_minus_less*100,color="#225b78")
    axes[1].scatter(arr.equal_pair_component_mean_more_minus_less*100,np.arange(len(arr)),color="#c16936",marker="D",label="Equal linked-group weight")
    axes[1].axvline(0,color="gray",lw=.8);axes[1].set(xlabel="More minus less of the feature: rank-spread change (points)",title="B. Arrangement associations are modest and context dependent")
    axes[1].legend(fontsize=8)
    for ext in ["png","pdf"]:fig.savefig(OUT/("sequence_order_and_arrangement."+ext),dpi=185,bbox_inches="tight")
    plt.close(fig)
    save(CHECKS,"sequence_analysis_checks")
    summary={"input_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"complete_sequences":len(z),"single_substitution_pairs":len(pairs),
             "substitution_components":int(pairs.component.nunique()),"composition_pair_count":len(order),
             "composition_linked_pair_components":int(order.pair_component.nunique()),"checks":len(CHECKS),
             "all_adjacent_swap_pairs_including_rotation_overlap":int(np.sum(all_adjacent)),
             "passed":sum(c["passed"] for c in CHECKS),"failed":[c for c in CHECKS if not c["passed"]],
             "elapsed_seconds":time.monotonic()-started,"definitions":{
                 "rank_range":"max minus min within-model percentile across the three families on the fixed full complete-case collection",
                 "centered_range":"log2 range after model-wise cohort median removal; prediction diagnostic, not calibrated MIC",
                 "local_reference":"exact expected score distance to another sequence matched on source and length, symmetric across anchors",
                 "component_subsets":"deterministic two-way split for sensitivity after exploration; not independent confirmation"},
             "python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,
             "dependency_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
                 ROOT/"scratchpad/source_sequences.tsv",ROOT/"analysis_review/outputs/sanity_checks/broader_sequence_components.csv",
                 ROOT/"analysis_review/outputs/additional_trends/exact_composition_two_member_pairs.csv"]},
             "scope":"Exploratory sequence/prediction associations; no measured activity, causal edit effects, or population significance inferred."}
    (OUT/"sequence_analysis_manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2));print(pd.DataFrame(opposite).to_string(index=False))
    assert all(c["passed"] for c in CHECKS)


if __name__=="__main__":main()
