#include "sltk_atom.h"
#include <iostream>

/*** Constructors and destructor ***/
FAtom::FAtom()
{
	d_x = 0.0;	
	d_y = 0.0;	
	d_z = 0.0;	
	type = 0;		
	natom = 0;
	cell_x = 0;
	cell_y = 0;
	cell_z = 0;
}
