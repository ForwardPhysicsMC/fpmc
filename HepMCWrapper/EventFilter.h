/*
 *  $Date: $
 *  $Revision: $
 *  
 */

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/PdfInfo.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>

class IsZ_Boson {
public:
  /// returns true if the GenParticle is a Z
  bool operator()( const HepMC::GenParticle* p ) { 
    if ( abs( p->pdg_id() ) == 23 ) return 1;
    return 0;
  }
};

class IsZZ_Event {

public:

   IsZZ_Event(bool debug = false): debug_(debug), idxZHad_(-1) { select_pid_ = false; }
   IsZZ_Event(std::vector<int> pid_list, bool debug = false): pid_list_( pid_list ), debug_(debug), idxZHad_(-1) { select_pid_ = false; }
   IsZZ_Event(int select_pid_1, int select_pid_2, bool debug = false): pid1_(select_pid_1), pid2_(select_pid_2), debug_(debug), idxZHad_(-1) { select_pid_ = true; }

   void SetHadronic( int idx ) { idxZHad_ = idx; }

   bool operator() ( const HepMC::GenEvent* evt ) {
       
      bool select = false;

      IsZ_Boson isZ;

      int n_Z = 0;
      std::vector<int> Z_decay_first_pid;
      std::vector<int> Z_decay_second_pid;
      for( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin() ; p != evt->particles_end() ; ++p ){

	 if ( isZ(*p) ){

	    if( (*p)->end_vertex() || (*p)->production_vertex() ) {
	       if( debug_ ) { std::cout << "======================" << std::endl; std::cout << "\t"; (*p)->print(); }
            }
	    if ( (*p)->end_vertex() ) {
	       if( debug_ ) { std::cout << "\t\t---Children" << std::endl; }
	       for ( HepMC::GenVertex::particle_iterator child = (*p)->end_vertex()->particles_begin(HepMC::children);
                                                         child != (*p)->end_vertex()->particles_end(HepMC::children); ++child ) {
		  if( debug_ ) { std::cout << "\t\t"; (*child)->print(); }
	       }
            }
	    if ( (*p)->production_vertex() ) {
	       if( debug_ ) { std::cout << "\t\t---Parents" << std::endl; }
	       for ( HepMC::GenVertex::particle_iterator parent = (*p)->production_vertex()->particles_begin(HepMC::parents);
                                                         parent != (*p)->production_vertex()->particles_end(HepMC::parents); ++parent ) {
		  if( debug_ ) { std::cout << "\t\t"; (*parent)->print(); }
	       }
	    }
	    if( (*p)->end_vertex() || (*p)->production_vertex() ) {
	       if( debug_ ) { std::cout << "======================" << std::endl; }
            }

            // Find decaying Z
	    if ( (*p)->end_vertex() && (*p)->production_vertex() ) {
               bool findParentZ = false; 
	       for ( HepMC::GenVertex::particle_iterator parent = (*p)->production_vertex()->particles_begin(HepMC::parents);
                                                         parent != (*p)->production_vertex()->particles_end(HepMC::parents); ++parent ) {
                  if( isZ(*parent) ) { findParentZ = true; break; }
	       }
               bool findChildZ = false; 
	       for ( HepMC::GenVertex::particle_iterator child = (*p)->end_vertex()->particles_begin(HepMC::children);
                                                         child != (*p)->end_vertex()->particles_end(HepMC::children); ++child ) {
                  if( isZ(*child) ) { findChildZ = true; break; }
	       }

               bool findLastZInChain = findParentZ && !findChildZ;
	       if( findLastZInChain ){
		  HepMC::GenVertex::particle_iterator it_part = (*p)->end_vertex()->particles_begin(HepMC::children);
		  HepMC::GenParticle* decay_part_1 = *it_part;
		  ++it_part;
		  HepMC::GenParticle* decay_part_2 = *it_part;

		  if( debug_ ){
		     std::cout << "Selected Z: "; (*p)->print();
		     std::cout << "           --> "; decay_part_1->print();       
		     std::cout << "           --> "; decay_part_2->print();       
		  }
          
		  int pdg_id_decay_1 = decay_part_1->pdg_id();
		  int pdg_id_decay_2 = decay_part_2->pdg_id();

		  Z_decay_first_pid.push_back( pdg_id_decay_1 );
		  Z_decay_second_pid.push_back( pdg_id_decay_2 );

                  ++n_Z;  
	       }
            }
	 }
      }

      // Select HepMC event
      if( select_pid_ ) {
	 if( n_Z == 2 ) {
	    if( ( abs( Z_decay_first_pid[0] ) == abs( pid1_ ) && 
		  abs( Z_decay_second_pid[0] ) == abs( pid1_ ) && 
		  Z_decay_first_pid[0]*Z_decay_second_pid[0] < 0 ) &&
		( abs( Z_decay_first_pid[1] ) == abs( pid2_ ) && 
		  abs( Z_decay_second_pid[1] ) == abs( pid2_ ) && 
		  Z_decay_first_pid[1]*Z_decay_second_pid[1] < 0 ) ) {
	       std::cout << "--- Selected HepMC event" << std::endl;
	       select = true;
	    }
	 }
      } else {
	 if( n_Z == 2 ) {
            bool select_Z1 = false;
            bool select_Z2 = false;
            if( idxZHad_ == 0 ) {
               if( abs( Z_decay_first_pid[0] ) >= 1 && abs( Z_decay_first_pid[0] ) <= 5 &&
                   Z_decay_second_pid[0] == -Z_decay_first_pid[0] ) select_Z1 = true;
            } else {
	       if( ( std::find( pid_list_.begin(), pid_list_.end(), abs( Z_decay_first_pid[0] ) ) != pid_list_.end() && 
		     Z_decay_second_pid[0] == -Z_decay_first_pid[0] ) ) select_Z1 = true;
            }

            if( idxZHad_ == 1 ) {
               if( abs( Z_decay_first_pid[1] ) >= 1 && abs( Z_decay_first_pid[1] ) <= 5 &&
                   Z_decay_second_pid[1] == -Z_decay_first_pid[1] ) select_Z2 = true;
            } else {
	       if( ( std::find( pid_list_.begin(), pid_list_.end(), abs( Z_decay_first_pid[1] ) ) != pid_list_.end() && 
	  	    Z_decay_second_pid[1] == -Z_decay_first_pid[1] ) ) select_Z2 = true; 
            }
            if( select_Z1 && select_Z2 ){
	       std::cout << "--- Selected HepMC event" << std::endl;
	       select = true;
            }
         }
      }

      return select;
   } // operator()

private:

   int pid1_;
   int pid2_;
   std::vector<int> pid_list_;

   bool select_pid_;
   bool debug_;
   int idxZHad_;

}; // class

