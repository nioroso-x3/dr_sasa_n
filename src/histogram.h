uint64
HST_GetBin(float,
           float,
           float,
           uint64,
           vector<float>);
class histogram
{
public:
  histogram();
  histogram(uint64,float,float); //bins, xmin, xmax
  histogram(const histogram& other);
  histogram& operator=(const histogram& other);
  void inc (float);
  void Pdist(void);
  void Initialize(uint64,float,float); //bins, xmin, xmax
  
  uint64 size;
  uint64 gcount;
  vector<float> range;
  vector<uint64> count;
  vector<float> p;
  ~histogram();
};
