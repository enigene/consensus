# Consensus v1.1, 12 Feb 2015
# Creating a consensus sequence from aligned set of FASTA sequences.
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
#
# v1.0, 11 Feb 2015 - first release
# v1.1, 12 Feb 2015 - added functions: arrMoveItemInPos, arrToPercenage, getConsensusStr
#
# Usage: gawk -v th=0.5 -v printmtx=1 -v printpmtx=1 -f cons.awk input-seq-aln.fas > output-cons.fas

BEGIN {
  DEL = "-";
  UNK = "N";

  # sets default threshold
  if ((th+0) <= 0) {
    th = sprintf("%.2f", 0.5);
  }
  if ((th+0) > 1) {
    th = sprintf("%.2f", th / 100);
  }

}

# makes two arrays with counted bases and dictionary
!/[^ACGTN-]/ {
  seq = $0;
  sLen = split(seq, seqA, "");
  for (i=1; i<=sLen; i++) {
    consA[i][seqA[i]]++;
    dictA[seqA[i]];
  }
}

END {
  dLen = asorti(dictA, dictAs);

  # move 'N' to the end of dict
  if (UNK in dictA) {
    arrMoveItemInPos(UNK, dLen, dictAs);
  }

  # move '-' to the end of dict
  if (DEL in dictA) {
    arrMoveItemInPos(DEL, dLen, dictAs);
  }

  # print consensus matrix
  if (printmtx) {
    # prints dict
    printf("#\t");
    for (i=1; i<=dLen; i++) {
      printf("%s", dictAs[i]);
      if (i < dLen) {
        printf("\t");
      }
    }
    printf("\n");
    # prints counted bases
    for (i=1; i<=sLen; i++) {
      printf("%d", i);
      for (j=1; j<=dLen; j++) {
        printf ("\t%d", consA[i][dictAs[j]]);
      }
      printf("\n");
    }
    printf("\n");
  }

  # makes percentage array
  arrToPercenage(consA);

  # prints percentage matrix
  if (printpmtx) {
    for (i=1; i<=sLen; i++) {
      printf("%d", i);
      for (j=1; j<=dLen; j++) {
        printf ("\t%.2f", consA[i][dictAs[j]]);
      }
      printf("\n");
    }
    printf("\n");
  }

  # finally gets consensus sequence
  consensus = getConsensusStr(consA, th, UNK);

  # prints consensus sequence
  printf("%s%.2f\n", ">cons_th_", th);
  print consensus;

}

# moves the named element into another position
function arrMoveItemInPos(value, toPos, array,    i, j, tmpArray, aLen) {
  i = j = 1;
  aLen = length(array);
  while (i <= aLen) {
    if (array[i] == value) {
      i++;
      continue;
    }
    if (j == toPos) {
      j++;
      continue;
    }
    tmpArray[j++] = array[i++];
  }
  tmpArray[toPos] = value;
  i = 1;
  while (i <= aLen) {
    array[i] = tmpArray[i];
    i++;
  }
}

# makes percentage array
function arrToPercenage(array,    i, j, aLen, sum) {
  aLen = length(array);
  for (i=1; i<=aLen; i++) {
    sum = 0;
    for (j in array[i]) {
      sum += array[i][j];
    }
    for (j in array[i]) {
      if (sum) {
        array[i][j] = sprintf("%.2f", array[i][j] / sum);
      }
    }
  }
}

function getConsensusStr(array, threshold, unkChar,    i, j, aLen, maxNum, currNum, consArr, consensus) {
  aLen = length(array);
  for (i=1; i<=aLen; i++) {
    maxNum = 0;
    # finds max base number in current pos
    for (j in array[i]) {
      currNum = array[i][j];
      if (currNum > maxNum) {
        maxNum = currNum;
      }
    }
    # makes new array with max base
    for (j in array[i]) {
      currNum = array[i][j];
      if (currNum == maxNum) {
        if (maxNum > sprintf("%.2f", threshold)) {
          consArr[i] = j;
        } else {
          consArr[i] = unkChar;
        }
      }
    }
    consensus = consensus consArr[i];
  }
  return consensus;
}

