#include "gtest/gtest.h"
#include <iomanip>
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "../relax.h"
#include "module_cell/unitcell.h"
#include "relax_test.h"
#include <fstream>


class Test_SETGRAD : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;
        Input_para& input = PARAM.input;
        UnitCell ucell;

        void SetUp()
        {
            PARAM.input.force_thr  = 0.001;
            PARAM.input.calculation = "cell-relax";

            ModuleBase::matrix force_in, stress_in;
            int nat = 3;

            force_in.create(nat,3);
            stress_in.create(3,3);

            force_in(0,0) = 0.1; force_in(0,1) = 0.1; force_in(0,2)= 0.1;
            force_in(1,0) = 0; force_in(1,1) = 0.1; force_in(1,2)= 0.1;
            force_in(2,0) = 0; force_in(2,1) = 0; force_in(2,2)= 0.1;

            stress_in(0,0) = 1; stress_in(0,1) = 1; stress_in(0,2)= 1;
            stress_in(1,0) = 0; stress_in(1,1) = 1; stress_in(1,2)= 1;
            stress_in(2,0) = 0; stress_in(2,1) = 0; stress_in(2,2)= 1;

            ucell.ntype = 1;
            ucell.nat = nat;
            ucell.atoms = new Atom[1];
            ucell.atoms[0].na = nat;
            ucell.omega = 1.0;
            ucell.lat0 = 1.0;
            
            ucell.iat2it = new int[nat];
            ucell.iat2ia = new int[nat];
            ucell.atoms[0].mbl.resize(nat);
            ucell.atoms[0].taud.resize(nat);
            ucell.atoms[0].tau.resize(nat);
            ucell.atoms[0].dis.resize(nat);
            ucell.lc = new int[3];

            ucell.iat2it[0] = 0;
            ucell.iat2it[1] = 0;
            ucell.iat2it[2] = 0;

            ucell.iat2ia[0] = 0;
            ucell.iat2ia[1] = 1;
            ucell.iat2ia[2] = 2;

            ucell.atoms[0].mbl[0].x = 0;
            ucell.atoms[0].mbl[0].y = 0;
            ucell.atoms[0].mbl[0].z = 1;

            ucell.atoms[0].mbl[1].x = 0;
            ucell.atoms[0].mbl[1].y = 1;
            ucell.atoms[0].mbl[1].z = 0;

            ucell.atoms[0].mbl[2].x = 1;
            ucell.atoms[0].mbl[2].y = 0;
            ucell.atoms[0].mbl[2].z = 0;

            ucell.atoms[0].taud[0] = 0.0;
            ucell.atoms[0].taud[1] = 0.0;
            ucell.atoms[0].taud[2] = 0.0;

            ucell.atoms[0].tau.resize(nat);

            ucell.lc[0] = 1;
            ucell.lc[1] = 1;
            ucell.lc[2] = 1;

            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);

            for(int i=0;i<3;i++)
            {
                result.push_back(ucell.atoms[0].taud[i].x);
                result.push_back(ucell.atoms[0].taud[i].y);
                result.push_back(ucell.atoms[0].taud[i].z);
            }
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            input.fixed_axes = "shape";
            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            input.fixed_axes = "volume";
            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);
            push_result();

            //reset lattice vector
            ucell.latvec.Identity();
            input.fixed_axes = "a"; //anything other than "None"
            input.fixed_ibrav = true;
            ucell.latName = "sc";
            ucell.lc[0] = 0;
            ucell.lc[1] = 0;
            ucell.lc[2] = 0;
            rl.init_relax(nat);
            rl.relax_step(ucell,force_in,stress_in,0.0);
            
            push_result();
        }

        void push_result()
        {
            result.push_back(ucell.latvec.e11);
            result.push_back(ucell.latvec.e12);
            result.push_back(ucell.latvec.e13);
            result.push_back(ucell.latvec.e21);
            result.push_back(ucell.latvec.e22);
            result.push_back(ucell.latvec.e23);
            result.push_back(ucell.latvec.e31);
            result.push_back(ucell.latvec.e32);
            result.push_back(ucell.latvec.e33);             
        }

};

TEST_F(Test_SETGRAD, relax_new)
{
    std::vector<double> result_ref = 
    {
        0,           0,            0.24293434145, 
        0,           0.242934341453,           0,
        0,           0,           0,
        //paramter for taud
        1.2267616333,0.2267616333,0.22676163333, 
        0,            1.2267616333 ,0.2267616333,
        0,            0,           1.22676163333,
        // paramter for fisrt time after relaxation
        1.3677603495, 0,            0,
        0,            1.36776034956, 0,
        0,            0,            1.36776034956,
        // paramter for second time after relaxation
        1.3677603495  ,0.3633367476,0.36333674766,
        0,            1.3677603495 ,0.36333674766,
        0,            0,            1.3677603495 ,
        // paramter for third time after relaxation
        1,0,0,0,1,0,0,0,1
        // paramter for fourth time after relaxation
    };
    for(int i=0;i<result.size();i++)
    {
        EXPECT_NEAR(result[i],result_ref[i],1e-8);
    }
}

class Test_RELAX : public testing::Test
{
    protected:
        Relax rl;
        std::vector<double> result;
        UnitCell ucell;

        void SetUp()
        {
            std::ifstream force_file("./support/force.txt");
            std::ifstream stress_file("./support/stress.txt");
            std::ifstream energy_file("./support/stress.txt");

            int nstep = 66;
            int nat = 5;
            double energy;
            PARAM.input.stress_thr = 0.01;

            this->setup_cell();

            ModuleBase::matrix force_in, stress_in;
            force_in.create(nat,3);
            stress_in.create(3,3);

            rl.init_relax(nat);

            for(int istep=0;istep<nstep;istep++)
            {
                for(int i=0;i<nat;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        force_file >> force_in(i,j);
                    }
                }
                for(int i=0;i<3;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        stress_file >> stress_in(i,j);
                    }
                }

                energy_file >> energy;

                PARAM.input.fixed_ibrav = false;
                rl.relax_step(ucell,force_in,stress_in,energy);

                result.push_back(ucell.atoms[0].taud[0].x);
                result.push_back(ucell.atoms[0].taud[0].y);
                result.push_back(ucell.atoms[0].taud[0].z);
                result.push_back(ucell.atoms[1].taud[0].x);
                result.push_back(ucell.atoms[1].taud[0].y);
                result.push_back(ucell.atoms[1].taud[0].z);
                result.push_back(ucell.atoms[2].taud[0].x);
                result.push_back(ucell.atoms[2].taud[0].y);
                result.push_back(ucell.atoms[2].taud[0].z);
                result.push_back(ucell.atoms[2].taud[1].x);
                result.push_back(ucell.atoms[2].taud[1].y);
                result.push_back(ucell.atoms[2].taud[1].z);
                result.push_back(ucell.atoms[2].taud[2].x);
                result.push_back(ucell.atoms[2].taud[2].y);
                result.push_back(ucell.atoms[2].taud[2].z);
                result.push_back(ucell.latvec.e11);
                result.push_back(ucell.latvec.e12);
                result.push_back(ucell.latvec.e13);
                result.push_back(ucell.latvec.e21);
                result.push_back(ucell.latvec.e22);
                result.push_back(ucell.latvec.e23);
                result.push_back(ucell.latvec.e31);
                result.push_back(ucell.latvec.e32);
                result.push_back(ucell.latvec.e33);
            }
        }

        void setup_cell()
        {
            int ntype = 3, nat = 5;
            ucell.ntype = ntype;
            ucell.nat = nat;

            ucell.omega = 452.590903143121;
            ucell.lat0 = 1.8897259886;
            ucell.iat2it = new int[nat];
            ucell.iat2ia = new int[nat];
            ucell.iat2it[0] = 0;
            ucell.iat2it[1] = 1;
            ucell.iat2it[2] = 2;
            ucell.iat2it[3] = 2;
            ucell.iat2it[4] = 2;

            ucell.iat2ia[0] = 0;
            ucell.iat2ia[1] = 0;
            ucell.iat2ia[2] = 0;
            ucell.iat2ia[3] = 1;
            ucell.iat2ia[4] = 2;

            ucell.atoms = new Atom[ntype];
            ucell.atoms[0].na = 1;
            ucell.atoms[1].na = 1;
            ucell.atoms[2].na = 3;
            
            for(int i=0;i<ntype;i++)
            {
                int na = ucell.atoms[i].na;
                ucell.atoms[i].mbl.resize(na);
                ucell.atoms[i].taud.resize(na);
                ucell.atoms[i].tau.resize(na);
                ucell.atoms[i].dis.resize(na);
                for (int j=0;j<na;j++)
                {
                    ucell.atoms[i].mbl[j] = {1,1,1};
                }
            }
            ucell.atoms[0].taud[0] = {0.5,0.5,0.00413599999956205};
            ucell.atoms[1].taud[0] = {0  ,0  ,0.524312999999893  };
            ucell.atoms[2].taud[0] = {0  ,0.5,0.479348999999274  };
            ucell.atoms[2].taud[1] = {0.5,0  ,0.479348999999274  };
            ucell.atoms[2].taud[2] = {0  ,0  ,0.958854000000429  };
            
            ucell.lc = new int[3];
            ucell.lc[0] = 1;
            ucell.lc[1] = 1;
            ucell.lc[2] = 1;

            ucell.latvec.e11 = 3.96;
            ucell.latvec.e12 = 0;
            ucell.latvec.e13 = 0;
            ucell.latvec.e21 = 0;
            ucell.latvec.e22 = 3.96;
            ucell.latvec.e23 = 0;
            ucell.latvec.e31 = 0;
            ucell.latvec.e32 = 0;
            ucell.latvec.e33 = 4.2768;
        }
};

TEST_F(Test_RELAX, relax_new)
{
    int size = 1584;
    double tmp;
    std::ifstream result_ref("./support/result_ref.txt");
    for(int i=0;i<size;i++)
    {
        result_ref >> tmp;
        EXPECT_NEAR(tmp,result[i],1e-8);
    }
}