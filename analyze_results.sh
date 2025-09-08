#!/usr/bin/env bash

RESULTS_DIR="results"
shopt -s nullglob
files=("$RESULTS_DIR"/*-run*.csv*)  # Note the * after .csv to handle double extensions
[[ ${#files[@]} -eq 0 ]] && { echo "No result files."; exit 1; }

printf "Dataset,BestValue,WorstValue,MeanValue,VarValue,MinTime,MaxTime,MeanTime,VarTime,Runs\n"

gawk -F',' '
  {
    # Extract dataset name, handling the double .csv issue
    split(FILENAME,p,"/");
    fname = p[length(p)];
    # Remove everything from -run onwards to get the base dataset name
    sub(/-run[0-9]+\.csv.*$/,"",fname);
    ds = fname;
    
    val=$2+0; t=$4+0            # +0 lets awk parse "8.1e+01"
    if (!(ds in best) || val<best[ds]) best[ds]=val;
    if (!(ds in worst) || val>worst[ds]) worst[ds]=val;
    sv[ds]+=val ; sv2[ds]+=val*val;
    st[ds]+=t   ; st2[ds]+=t*t;
    if (!(ds in mint) || t<mint[ds]) mint[ds]=t;
    if (!(ds in maxt) || t>maxt[ds]) maxt[ds]=t;
    n[ds]++;
  }
  END{
    PROCINFO["sorted_in"]="@ind_str_asc";
    for (ds in n){
      mv=sv[ds]/n[ds];  mt=st[ds]/n[ds];
      vv=sv2[ds]/n[ds]-mv*mv;
      vt=st2[ds]/n[ds]-mt*mt;
      printf "%s,%.10e,%.10e,%.10e,%.10e,%.4f,%.4f,%.4f,%.4f,%d\n",
             ds,best[ds],worst[ds],mv,vv,mint[ds],maxt[ds],mt,vt,n[ds];
    }
  }' "${files[@]}"
