struct NUCL{
  string nuc;
  NUCL *next;
};

struct PROPERTY;

struct GENELC{
  string g;
  string odir;
  int ex;
  int cnt;
  int npcnt;
  int egcnt;
  int mglen;
  int mnlen;
  int *melen;
  GENELC *next;
  ofstream outl,outg,outgn,outn,outu5,outu3,outall,outgnutr;
  ofstream outog,outogn;
  ofstream *oute,*outi;
  ofstream *outien,*outin,*outienl,*outienr; //28-Aug-2014
  ofstream *outen;
  PROPERTY **sameex,**samein,**sameutr;
};

struct PROPERTY{
  string g;
  string name;
  string al;
  bool par;
  bool pse,unuse;
  bool eg;
  int a[4];  
  int excnt;
  int length;
  int *es,*ee;
  int *eid;
  string trans;
  NUCL *ntop;
  PROPERTY *next;
  char *seq;
  GENELC *glc;
};


