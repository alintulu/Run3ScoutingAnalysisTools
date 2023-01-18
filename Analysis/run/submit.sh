#!/bin/bash

while read LINE; do
  DATASET=$LINE
  NAME=$(echo "$LINE" | cut -d "/" -f2)
  PROCESS=$(echo "$NAME" | cut -d "_" -f1)
  if [[ $LINE == *"Run3Summer22EE"* ]]; then
    CAMPAIGN="Run3Summer22EE"
  else
    CAMPAIGN="Run3Summer22"
  fi

  if [[ -f "${PROCESS}/${CAMPAIGN}/${NAME}.py" ]]; then
    echo "Already submitted ${PROCESS}/${CAMPAIGN}/${NAME}.py"
    continue
  fi

  echo DATASET=$DATASET
  echo NAME=$NAME
  echo PROCESS=$PROCESS
  echo CAMPAIGN=$CAMPAIGN

  mkdir -p ${PROCESS}/${CAMPAIGN}
  cp template.py ${PROCESS}/${CAMPAIGN}/${NAME}.py
  sed -i "s@DATASET@$DATASET@" ${PROCESS}/${CAMPAIGN}/${NAME}.py
  sed -i "s@NAME@$NAME@" ${PROCESS}/${CAMPAIGN}/${NAME}.py
  sed -i "s@CAMPAIGN@$CAMPAIGN@" ${PROCESS}/${CAMPAIGN}/${NAME}.py

  crab submit ${PROCESS}/${CAMPAIGN}/${NAME}.py 
done < files.txt
