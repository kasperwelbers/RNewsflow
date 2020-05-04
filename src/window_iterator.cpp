#include "window_iterator.h"

void windowIterator::start(const SpMat& m) {
  if (!(m.cols() == index2.size())) stop("Number of columns in m needs to match the size of index2");
  rightsum = std::vector<double>(m.rows());
  leftsum = std::vector<double>(m.rows());
  int date, first_date, loop_date;
  
  // find first position that has a full left window
  first_date = std::get<1>(index2[0]);
  for (; pos < index1.size(); pos++) {
    date = std::get<1>(index1[pos]);
    if (first_date - date < lwindow) break;
  }
  
  // find the ls, le, rs and re positions by looping over index2 and comparing date to pos date
  bool ls_flag = false, rs_flag = false, done_flag = false;
  for (int j = 0; j < index2.size(); j++) {
    loop_date = std::get<1>(index2[j]);
    // first time in left window is ls
    if (loop_date - date > lwindow) {
      if (!ls_flag) {
        ls = j;
        ls_flag = true;
      }
    }
    // until le is reached, sum up the left side vocabulary
    if (loop_date < (date + lwindow_border)) {
      le = j;
      for (SpMat::InnerIterator itleft(m, j); itleft; ++itleft) {
        leftsum[itleft.row()] += itleft.value();
        lefttotal += itleft.value();
      }
    } 
    // once rs is reached, sum up the right side vocabulary 
    if (loop_date > (date + rwindow_border)) {
      if (!rs_flag) {
        rs = j;
        rs_flag = true;
      }
      for (SpMat::InnerIterator itright(m, j); itright; ++itright) {
        rightsum[itright.row()] += itright.value();
        righttotal += itright.value();
      }
    }  
    // stop once re is reached
    if (loop_date - date > rwindow) {
      re = j - 1;
      done_flag = true;
      break;
    }
  }
  if (!done_flag) stop("date range is too small for at least one full window (left and right) to occur");
}

void windowIterator::increment(const SpMat& m) {
  if (pos == index1.size() - 1) {
    done=true;
    return;
  }
  
  pos++;
  int date = std::get<1>(index1[pos]);
  for (;ls < index2.size(); ls++) {
    int lsdate = std::get<1>(index2[ls]);
    if (lsdate - date > lwindow) break;
    for (SpMat::InnerIterator itleft(m, ls); itleft; ++itleft) {
      leftsum[itleft.row()] -= itleft.value();
      lefttotal -= itleft.value();
    }
  }
  for (;le < index2.size(); le++) {
    int ledate = std::get<1>(index2[le]);
    if (ledate >= (date + lwindow_border)) break;
    for (SpMat::InnerIterator itleft(m, le); itleft; ++itleft) {
      leftsum[itleft.row()] += itleft.value();
      lefttotal += itleft.value();
    }
  }
  for (;rs < index2.size(); rs++) {
    int rsdate = std::get<1>(index2[rs]);
    if (rsdate >= (date + rwindow_border)) break;
    for (SpMat::InnerIterator itright(m, rs); itright; ++itright) {
      rightsum[itright.row()] -= itright.value();
      righttotal -= itright.value();
    }
  }
  for (;re < index2.size(); re++) {
    int redate = std::get<1>(index2[re]);
    if (redate - date > rwindow) break;
    if (re == index2.size() - 1) done=true; // if end of right window (re) is not in sight (beyond last value in index2), we're done
    for (SpMat::InnerIterator itright(m, re); itright; ++itright) {
      rightsum[itright.row()] += itright.value();
      righttotal += itright.value();
    }
  }
}