#include "read_stru.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"

namespace unitcell
{
    bool check_tau(const Atom* atoms,
                   const int ntype) 
    {
        ModuleBase::TITLE("UnitCell","check_tau");
        ModuleBase::timer::tick("UnitCell","check_tau");
        
        ModuleBase::Vector3<double> diff = 0.0;
        double norm = 0.0;
        double tolerence_bohr = 1.0e-3;

        //GlobalV::ofs_running << "\n Output nearest atom not considering periodic boundary condition" << std::endl;
        //GlobalV::ofs_running << " " << std::setw(5) << "TYPE" << std::setw(6) << "INDEX" 
        //<< std::setw(20) << "NEAREST(Bohr)" 
        //<< std::setw(20) << "NEAREST(Angstrom)" << std::endl; 
        for(int T1=0; T1< ntype; T1++)
        {
            for(int I1=0; I1< atoms[T1].na; I1++)
            {    
                double shortest_norm = 10000.0; // a large number
                //int nearest_atom_type = 0;
                //int nearest_atom_index = 0;
                for(int T2=0; T2<ntype; T2++)
                {
                    for(int I2=0; I2<atoms[T2].na; I2++)
                    {
                        if(T1==T2 && I1==I2)
                        {
                            shortest_norm = 0.0;
                            //nearest_atom_type = T1;
                            //nearest_atom_index = I2;
                            // self atom
                        }
                        else
                        {
                            diff = atoms[T1].tau[I1] - atoms[T2].tau[I2];
                            norm = diff.norm() * lat0;
                            if( shortest_norm > norm )
                            {
                                shortest_norm = norm;
                                //nearest_atom_type = T2;
                                //nearest_atom_index = I2;
                            }
                            if( norm < tolerence_bohr ) // unit is Bohr
                            {    
                                GlobalV::ofs_warning << " two atoms are too close!" << std::endl;
                                GlobalV::ofs_warning << " type:" << atoms[T1].label << " atom " << I1 + 1 << std::endl; 
                                GlobalV::ofs_warning << " type:" << atoms[T2].label << " atom " << I2 + 1 << std::endl; 
                                GlobalV::ofs_warning << " distance = " << norm << " Bohr" << std::endl;
                                return false;
                            }
                        }
                    }
                }
                //GlobalV::ofs_running << " " << std::setw(5) << atoms[T1].label << std::setw(6) << I1+1 
                //<< std::setw(20) << shortest_norm  
                //<< std::setw(20) << shortest_norm * ModuleBase::BOHR_TO_A << std::endl;
            }
        }

        ModuleBase::timer::tick("UnitCell","check_tau");
        return true;
    }
}