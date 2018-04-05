#include "stdafx.h"
#include "histogram.h"

uint64
HST_GetBin(float xmin,
           float xmax,
           float value,
           uint64 bins,
           vector<float> binrange){
  if(bins == 0) return 0;
  if(value >= binrange[bins-1]) return bins - 1;
  if(value < xmin) return 0;
  for (uint64 i = 0; i < bins; ++i){
    if((value >= binrange[i]) && (value < binrange[i+1]))
      return i;
  }
  throw runtime_error("HST_GetBin: Invalid data?\n");
  
}

histogram::histogram(){
}

histogram::histogram(uint64 bins,
                     float xmin,
                     float xmax){
  float step = (xmax - xmin) / bins;
  gcount = 0;
  size = bins;
  count.resize(size,0);
  p.resize(size,0);
  range.resize(size+1,0);
  for (uint64 i = 0; i < (size+1); i++){
    range[i] = (xmin+i*step);
  }
}

histogram&
histogram::operator=(const histogram& other){
  size = other.size;
  gcount = other.gcount;
  range = other.range;
  count = other.count;
  p = other.p;
  return *this;
}

histogram::histogram(const histogram& other){
  size = other.size;
  gcount = other.gcount;
  range = other.range;
  count = other.count;
  p = other.p;
}

void
histogram::Initialize(uint64 bins,
                      float xmin,
                      float xmax){
  float step = (xmax - xmin) / bins;
  gcount = 0;
  size = bins;
  count.resize(size,0);
  p.resize(size,0);
  range.resize(size+1,0);
  for (uint64 i = 0; i < (size+1); i++){
    range[i] = (xmin+i*step);
  }
}


void
histogram::inc (float value){
  //cout << value << " -> ";
  if (range.size() == 0) throw runtime_error("histogram::inc Initialize histogram first");
  else{
    if( value >= range.back()){
        gcount++;
        count.back()++;
    }
    else{
      for (uint64 i = 0; i < range.size() - 1; ++i){
        if(value >= range[i] && value < range[i+1]){
          gcount++;
          count[i]++;
          break;
        }
      }
    }
  }
}

void
histogram::Pdist(void){
  if(gcount > 0){
    for (uint64 i = 0; i < count.size(); i++){
      p[i] = (float)count[i] / (float)gcount;
    }
  }
  else{ 
    cerr << "Empty histogram!\n";
  }
}

histogram::~histogram()
{
}
