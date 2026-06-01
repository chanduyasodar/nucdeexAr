#ifndef __NUCDEEXNUCLEUSTABLE__HH__
#define __NUCDEEXNUCLEUSTABLE__HH__

#include <map>

#include "NucDeExNucleus.hh"

class NucDeExNucleusTable{
  public:
  NucDeExNucleusTable();
  virtual ~NucDeExNucleusTable(){;};

  bool ReadTables(const bool init_flag=1);
    // 1 -> new arraies
  int getID(const char* name);
  int GetNumofNuc(){return num_of_nuc;};
  
  NucDeExNucleus* GetNucleusPtr(const char* name);
  NucDeExNucleus* GetNucleusPtr(int id);
  NucDeExNucleus* GetNucleusPtr(int Z,int N);
  NucDeExNucleus* GetNucleusPtrPDG(int PDG);
  const char* nuc_name[119]
    = {"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
 "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
 "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
 "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
 "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
 "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
 "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
 "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
 "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
 "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"}
; // [Z]

  private:
  int num_of_nuc;
  bool flag_read;
  NucDeExNucleus* _nucleus;
  std::map<std::string, int> _nucleus_id;
  std::map<std::string, int> :: iterator _p_id;
};
#endif
