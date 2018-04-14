
#include "parameter.cuh"

parameter::parameter(char* inputdir):rf()
{
	input_file_name = new char [1024];

	strcpy(input_file_name, inputdir);
    rf.openinput( inputdir ); 

    initial_particle();
	initial_grid();
	initial_static_potential();
	inital_time_set();
};

//////////////////////////////////////////////////////////////////////////////////////////

void parameter::initial_particle( void )
// read and particle information
{

  rf.openinput( input_file_name );

      energy_e      = atof( rf.setget( "&particle", "energy_e" ) );
      energy_ion    = atof( rf.setget( "&particle", "energy_ion"  ) );
      Ib            = atof( rf.setget( "&particle", "Ib" ) );
      wegith        = atof( rf.setget( "&particle", "weigth" ) );
      ni            = atof( rf.setget( "&particle", "atom_density" ) );

  rf.closeinput();

}

//////////////////////////////////////////////////////////////////////////////////////////
void parameter::initial_grid( void )
// read and gridinformation
{

  rf.openinput( input_file_name );

     dr        = atof( rf.setget( "&box", "cells_dr" ) );
     dz        = atof( rf.setget( "&box", "cells_dz"  ) );
	 L         = atof( rf.setget( "&box", "L" ) );
     R         = atof( rf.setget( "&box", "R" ) );	
     ca_rd     = atof( rf.setget( "&box", "cathode_rad" ) );
	 ca_len    = atof( rf.setget( "&box", "cathode_len" ) );
     rf.closeinput();

}
void parameter::initial_static_potential( void )
// read and cell information
{

  rf.openinput( input_file_name );

     anode_p      = atof( rf.setget( "&static_field", "wall_electric_potential" ) );
     cathode_p    = atof( rf.setget( "&static_field", "cathode_electric_potential"  ) );
     screen_p     = atof( rf.setget( "&static_field", "screen_electric_potential" ) );

  rf.closeinput();

}

void parameter::inital_time_set( void )
// read and cell information
{

  rf.openinput( input_file_name );

     total_time      = atof( rf.setget( "&time_set", "total_time" ) );
     time_inter      = atof( rf.setget( "&time_set", "time_interval"  ) );

  rf.closeinput();

}


//EOF

