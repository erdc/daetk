#include "ConstBCECDM.h"

namespace Daetk 
{

using std::cerr;
using std::cout;
using std::endl;
using std::remove_if;

void ConstBCECDM::setBoundaryConditions(face_type ft, 
					condition_type ct, 
					Petsc::StencilMM& node,
					real bc_value, 
					Vec* var, 
					Vec* varprime,
					real* value1,
					real* value2) 
{
  //should be able just to use isLocal here since 
  //dirichle values should be loaded into solution with ECDM
  if (node.isInGhostRegion)
    {
      switch (ct)
        {
        case DIRICHLET:
          {
            assert(var != 0);
            assert(varprime != 0);
            Dirichlet d={ft,node.petscToGlobal(node.center),bc_value,
			 value1,var,varprime};
            dirichlet.push_back(d);
            break;
          }
        case NEUMANN:
          {
            Neumann n={ft,node.petscToGlobal(node.center),bc_value};
            neumann.push_back(n);
            break;
          }
        case ROBBINS:
          {
            assert(var != 0);
            assert(varprime != 0);
	    //value1 goes to dirMult, value2 goes to neuMult
            Robbins r={ft,node.petscToGlobal(node.center),bc_value,
			 *value1,*value2,var,varprime};
            robbins.push_back(r);
            break;
          }
        case INTERIOR:
          {
            InteriorBC in={ft,node.petscToGlobal(node.center)};
            interior.push_back(in);
            break;
          }
        case DUMMY:
          {
            assert(var != 0);
            assert(varprime != 0);
            DummyBC du={ft,node.petscToGlobal(node.center),bc_value,
			 var,varprime};
            dummy.push_back(du);
            break;
          }
        case GRADIENT:
          {
            Gradient g={ft,node.petscToGlobal(node.center),bc_value};
            noflow.push_back(g);
            break;
          }
        default:
          {
            cerr<<"Boundary condition type "<<ft<<"not implemented"<<endl;
            exit(1);
          }
        }
    }
}

void ConstBCECDM::print()
{
  //try and make sure no duplicates here for now
  checkForDuplicateBCs();

  dit=dirichlet.begin();
  while(dit!=dirichlet.end())
    {
      cout<<"D"<<'\t'<<dit->face<<'\t'<<dit->n<<'\t'<<dit->value;
      if (dit->scale)
        cout<<'\t'<<(*dit->scale)<<endl;
      else
        cout<<endl;
      ++dit;
    }
  
  nit=neumann.begin();
  while(nit!=neumann.end())
    {
      cout<<"N"<<'\t'<<nit->face<<'\t'<<nit->n<<'\t'<<nit->value<<endl;
      ++nit;
    }

  ruit=robbins.begin();
  while(ruit!=robbins.end())
    {
      cout<<"R"<<'\t'<<ruit->face<<'\t'<<ruit->n<<'\t'<<ruit->value;
      if (ruit->dirMult && ruit->neuMult)
        cout<<'\t'<<(ruit->dirMult)<<'\t'<<(ruit->neuMult)<<endl;
      else
        cout<<endl;
      ++ruit;
    }

  init=interior.begin();
  while(init!=interior.end())
    {
      cout<<"IN"<<'\t'<<init->face<<'\t'<<init->n<<'\t'<<endl;
      ++init;
    }

  duit=dummy.begin();
  while(duit!=dummy.end())
    {
      cout<<"DU"<<'\t'<<duit->face<<'\t'<<duit->n<<'\t'<<duit->value
	  <<'\t'<<endl;
      ++duit;
    }

  noit=noflow.begin();
  while(noit!=noflow.end())
    {
      cout<<"G"<<'\t'<<noit->face<<'\t'<<noit->n<<'\t'
	  <<noit->value<<endl;
      ++noit;
    }
}

//mwf try and make sure that no duplication between dummy list and others
void ConstBCECDM::checkForDuplicateBCs()
{
  duit=dummy.begin();
  while(duit!=dummy.end())
    {
      int dummyNode = duit->n;

      //check Neumann first
      BCnodeIsTheSame<Neumann> neuComp(dummyNode);
      std::vector<Neumann>::iterator newNeuEnd =
	remove_if(neumann.begin(),neumann.end(),neuComp);

      neumann.erase(newNeuEnd,neumann.end());

      //next try Dirichlet
      BCnodeIsTheSame<Dirichlet> dirComp(dummyNode);
      std::vector<Dirichlet>::iterator newDirEnd =
	remove_if(dirichlet.begin(),dirichlet.end(),dirComp);

      dirichlet.erase(newDirEnd,dirichlet.end());

      //next try Robbins
      BCnodeIsTheSame<Robbins> robComp(dummyNode);
      std::vector<Robbins>::iterator newRobEnd =
	remove_if(robbins.begin(),robbins.end(),robComp);

      robbins.erase(newRobEnd,robbins.end());

      //next try Gradient
      BCnodeIsTheSame<Gradient> noComp(dummyNode);
      std::vector<Gradient>::iterator newNoEnd =
	remove_if(noflow.begin(),noflow.end(),noComp);

      noflow.erase(newNoEnd,noflow.end());

      //last check interior
      BCnodeIsTheSame<InteriorBC> intComp(dummyNode);
      std::vector<InteriorBC>::iterator newIntEnd =
	remove_if(interior.begin(),interior.end(),intComp);

      interior.erase(newIntEnd,interior.end());

      
      ++duit;
    }//end duit while
}//end function

//try and do alot of common operations here now
//========== 1d ==========
void 
ConstBCECDM::applyBoundaryConditionsToPressureDifferences(Petsc::StencilMM& 
							  node,
							  Vec& pDiff_x)
{

  //==============================dirichlet ======================
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_x[node.interLeftGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==============================robbins======================
  
  ruit = robbins.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"R"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop

  //==================== neumann======================
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      //mwf make this isInGhostRegion?
      if (node.isInGhostRegion)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
 
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		pDiff_x[node.interLeftGhost] = 0.0;
		//              cout<<"FL"<<'\t'<<flux_x[node.interLeft]<<'\t'<<nit->value<<'\t'<< flux_x[node.interRight]<<endl;
		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		//              cout<<"FR"<<'\t'<<flux_x[node.interRight] <<'\t'<<nit->value <<'\t'<<flux_x[node.interLeft]<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }//end nit loop

  //====================dummy===========================

  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(duit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    } //end duit loop
  
  //====================interior===========================

  //==================== no flow======================
  //still just zero these for now, may want to try and reuse
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //mwf make this isInGhostRegion?
      if (node.isInGhostRegion)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
 
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		pDiff_x[node.interLeftGhost] = 0.0;
		//              cout<<"FL"<<'\t'<<flux_x[node.interLeft]<<'\t'<<nit->value<<'\t'<< flux_x[node.interRight]<<endl;
		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		//              cout<<"FR"<<'\t'<<flux_x[node.interRight] <<'\t'<<nit->value <<'\t'<<flux_x[node.interLeft]<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end noit loop

}//end function


//========== 2d ==========
void 
ConstBCECDM::applyBoundaryConditionsToPressureDifferences(Petsc::StencilMM& 
							  node,
							  Vec& pDiff_x,
							  Vec& pDiff_y)
{

  //==============================dirichlet ======================
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_x[node.interLeftGhost] = 0.0;
		pDiff_y[node.interBackGhost]  = 0.0;
		pDiff_y[node.interFrontGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost] = 0.0;
		pDiff_y[node.interBackGhost]  = 0.0;
		pDiff_y[node.interFrontGhost] = 0.0;
		break;
	      }
	    case FRONT:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_y[node.interFrontGhost] = 0.0;
		pDiff_x[node.interRightGhost] = 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    case BACK:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==============================robbins======================
  
  ruit = robbins.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_x[node.interLeftGhost] = 0.0;
		pDiff_y[node.interBackGhost]  = 0.0;
		pDiff_y[node.interFrontGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost] = 0.0;
		pDiff_y[node.interBackGhost]  = 0.0;
		pDiff_y[node.interFrontGhost] = 0.0;
		break;
	      }
	    case FRONT:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_y[node.interFrontGhost] = 0.0;
		pDiff_x[node.interRightGhost] = 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    case BACK:
	      {
		//mwf should gravity term be different for boundary?
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"R"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop

  //==================== neumann======================
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      //mwf make this isInGhostRegion?
      if (node.isInGhostRegion)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
 
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		pDiff_x[node.interLeftGhost] = 0.0;
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_y[node.interFrontGhost]= 0.0;
		
		//              cout<<"FL"<<'\t'<<flux_x[node.interLeft]<<'\t'<<nit->value<<'\t'<< flux_x[node.interRight]<<endl;
		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_y[node.interFrontGhost]= 0.0;
		//              cout<<"FR"<<'\t'<<flux_x[node.interRight] <<'\t'<<nit->value <<'\t'<<flux_x[node.interLeft]<<endl;
		break;
	      }
	    case FRONT:
	      {
		pDiff_y[node.interFrontGhost]= 0.0;
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;

		//              cout<<"FF"<<'\t'<<flux_y[node.interFront]<<'\t'<<nit->value <<'\t'<< flux_y[node.interBack]<<endl;
		break;
	      }
	    case BACK:
	      {
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		//              cout<<"FB"<<'\t'<<flux_y[node.interBack] <<'\t'<<nit->value<<'\t'<< flux_y[node.interFront]<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }//end nit loop

  //====================dummy===========================

  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(duit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_y[node.interFrontGhost]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    } //end duit loop
  
  //====================interior===========================

  init = interior.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_y[node.interFrontGhost]= 0.0;
		break;
	      }
	    case INTERIOR_Y:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<init->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++init;
    } //end init loop

  //==================== no flow======================
  //still just zero these for now, may want to try and reuse
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //mwf make this isInGhostRegion?
      if (node.isInGhostRegion)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
 
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		pDiff_x[node.interLeftGhost] = 0.0;
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_y[node.interFrontGhost]= 0.0;
		
		//              cout<<"FL"<<'\t'<<flux_x[node.interLeft]<<'\t'<<nit->value<<'\t'<< flux_x[node.interRight]<<endl;
		break;
	      }
	    case RIGHT:
	      {
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_y[node.interFrontGhost]= 0.0;
		//              cout<<"FR"<<'\t'<<flux_x[node.interRight] <<'\t'<<nit->value <<'\t'<<flux_x[node.interLeft]<<endl;
		break;
	      }
	    case FRONT:
	      {
		pDiff_y[node.interFrontGhost]= 0.0;
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;

		//              cout<<"FF"<<'\t'<<flux_y[node.interFront]<<'\t'<<nit->value <<'\t'<< flux_y[node.interBack]<<endl;
		break;
	      }
	    case BACK:
	      {
		pDiff_y[node.interBackGhost] = 0.0;
		pDiff_x[node.interRightGhost]= 0.0;
		pDiff_x[node.interLeftGhost] = 0.0;
		//              cout<<"FB"<<'\t'<<flux_y[node.interBack] <<'\t'<<nit->value<<'\t'<< flux_y[node.interFront]<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end noit loop

}//end function

//this definitely won't work for multiple variables
//I don't think I need this one now
//================== 1d =======================
void 
ConstBCECDM::applyBoundaryConditionsToFluxes(Petsc::StencilMM& 
					     node,
					     Vec& flux_x)
{

  //==================== dirichlet ====================
  
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==================== robbins ====================
  
  ruit = robbins.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop

  //==================== neumann =======================

  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      //mwf now have to add isLocal because I'm keeping bc's for 
      //mwf ghost region
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
 
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }//end nit loop

  //==================== dummy ======================
  
  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(duit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    } //end duit loop


  //==================== no flow =======================
  //mwf just zero these for now, could try and load correct value into
  //mwf flux at exterior?
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //mwf now have to add isLocal because I'm keeping bc's for 
      //mwf ghost region
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
 
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end noit loop

}//end function


//this definitely won't work for multiple variables
//I don't think I need this one now
//================== 2d =======================
void 
ConstBCECDM::applyBoundaryConditionsToFluxes(Petsc::StencilMM& 
					     node,
					     Vec& flux_x,
					     Vec& flux_y)
							
{

  //==================== dirichlet ====================
  
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case FRONT:
	      {
		flux_y[node.interFront]= 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case BACK:
	      {
		flux_y[node.interBack] = 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==================== robbins ====================
  
  ruit = robbins.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case FRONT:
	      {
		flux_y[node.interFront]= 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case BACK:
	      {
		flux_y[node.interBack] = 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop

  //==================== neumann =======================
  //mwf next try and set calculated flux into "outside flux"
  //mwf and then set inside flux to actual value?
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      //mwf now have to add isLocal because I'm keeping bc's for 
      //mwf ghost region
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
 
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		//interior node is calculated correctly for outside flux
		//flux_x[node.interLeft] = flux_x[node.interRight];
		//flux_x[node.interRight]= -nit->value;
		flux_x[node.interLeft] = 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case RIGHT:
	      {
		//flux_x[node.interRight]= flux_x[node.interLeft];
		//flux_x[node.interLeft] = nit->value;
		flux_x[node.interRight]= 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case FRONT:
	      {
		//flux_y[node.interFront]= flux_y[node.interBack];
		//flux_y[node.interBack] = -nit->value;
		flux_y[node.interFront]= 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case BACK:
	      {
		//flux_y[node.interBack] = flux_y[node.interFront];
		//flux_y[node.interFront]= nit->value;
		flux_y[node.interBack] = 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }//end nit loop

  //==================== dummy ======================
  
  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(duit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    } //end duit loop

  //==================== interior ========================
  
  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		//              cout<<"FL"<<'\t'<<flux_x[node.interLeft]<<'\t'<<nit->value<<'\t'<< flux_x[node.interRight]<<endl;
		break;
	      }
	    case INTERIOR_Y:
	      {
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		//              cout<<"FL"<<'\t'<<flux_x[node.interLeft]<<'\t'<<nit->value<<'\t'<< flux_x[node.interRight]<<endl;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h:Face Type "<<init->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++init;
    } //end duit loop

  //==================== no flow =======================
  //mwf just zero these for now, could try and load correct value into
  //mwf flux at exterior?
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //mwf now have to add isLocal because I'm keeping bc's for 
      //mwf ghost region
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
 
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		//where should I account for the sign of normal vector??
		flux_x[node.interLeft] = 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case RIGHT:
	      {
		flux_x[node.interRight]= 0.0;
		flux_y[node.interBack] = 0.0;
		flux_y[node.interFront]= 0.0;
		break;
	      }
	    case FRONT:
	      {
		flux_y[node.interFront]= 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    case BACK:
	      {
		flux_y[node.interBack] = 0.0;
		flux_x[node.interRight]= 0.0;
		flux_x[node.interLeft] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end noit loop

}//end function

//================ 1d ============
void 
ConstBCECDM::applyBoundaryConditionsToPressureDerivatives(Petsc::StencilMM& 
							  node,
							  Vec& DpDiff_x_center,
							  Vec& DpDiff_x_right)
{
  //==============================dirichlet======================

  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==================== robbins ======================

  ruit = robbins.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"R"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop

  //==================================neumann======================
  nit = neumann.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++nit;
    } //end nit loop
  //==============================dummy======================

  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(duit->n);
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {

		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    } //end duit loop

  //==================================no flow======================
  //mwf could try and shift these over too
  noit = noflow.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++noit;
    } //end noit loop

}//end function



//only works for 1 variable
//================ 2d ============
void 
ConstBCECDM::applyBoundaryConditionsToPressureDerivatives(Petsc::StencilMM& 
							  node,
							  Vec& DpDiff_x_center,
							  Vec& DpDiff_x_right,
							  Vec& DpDiff_y_center,
							  Vec& DpDiff_y_back)
{
  //==============================dirichlet======================

  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case FRONT:
	      {
		DpDiff_y_center[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    case BACK:
	      {
		DpDiff_y_back[node.interBackGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==================== robbins ======================

  ruit = robbins.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case FRONT:
	      {
		DpDiff_y_center[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    case BACK:
	      {
		DpDiff_y_back[node.interBackGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"R"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop

  //==================================neumann======================
  nit = neumann.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case FRONT:
	      {
		DpDiff_y_center[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    case BACK:
	      {
		DpDiff_y_back[node.interBackGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++nit;
    } //end nit loop
  //==============================dummy======================

  duit = dummy.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(duit->n);
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {

		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;
		DpDiff_y_center[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++duit;
    } //end duit loop

  //========================interior======================
  init = interior.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		DpDiff_y_center[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interBackGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;


		break;
	      }
	    case INTERIOR_Y:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interLeftGhost]= 0.0;

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<init->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
      //cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++init;
    } //end init loop
  //==================================no flow======================
  //mwf could try and shift these over too
  noit = noflow.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      //mwf only set if node is local to this processor?
      //      if (node.isLocal)
      //mwf now try to set bc's for ghost region too so that I can 
      //mwf get all the transverse derivatives
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		DpDiff_x_center[node.interLeftGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		DpDiff_x_right[node.interRightGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;

		break;
	      }
	    case FRONT:
	      {
		DpDiff_y_center[node.interFrontGhost] = 0.0;
		DpDiff_y_back[node.interFrontGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    case BACK:
	      {
		DpDiff_y_back[node.interBackGhost] = 0.0;
		DpDiff_y_center[node.interBackGhost] = 0.0;
		DpDiff_x_right[node.interLeftGhost] = 0.0;
		DpDiff_x_center[node.interRightGhost]= 0.0;

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++noit;
    } //end noit loop

}//end function


//============== 1d ======================
//this definitely won't work for multiple variables
//haven't looked at this one yet
void ConstBCECDM::applyBoundaryConditionDerivativesECDM(Petsc::StencilMM& node,
							Vec& Dflux_x_center, 
							Vec& Dflux_x_right)
{
  
  //====================dirichlet===========================
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==================== robbins ===========================
  ruit = robbins.begin();
  //std::cout<<ruit<<'\t'<<robbins.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop


  //====================neumann===========================
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      if (node.isLocal)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }//end nit loop

  //====================== dummy =========================
  //zero everything
  duit = dummy.begin();
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      if (node.isLocal)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(duit->n);
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++duit;
    } //end duit loop

  //==================== interior ====================
  //zero transverse
  //and add in extra term for interface coefficient dependence
  //on left/right front/back neighbor
  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf now have to add isLocal because I'm keeping bc's for 
      //mwf ghost region
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
#ifdef USE_CONT_INTERFACE_TERMS_ECDM
		assert(currentDiv);

		//now get extra term
		real fluxDeriv(0.0);
		int face(0);
		//interface coefficient is upwind(Kr_{left},Kr_{right})*
		//1/2(rho_{left}+rho_{right})
	      
		//flux is from the right and so gets a negative
		//even though derivative is wrt left variable
		face = 1;
		fluxDeriv = currentDiv->getInteriorFluxDerivTerm_x(face);
		Dflux_x_center[node.interLeft] += -fluxDeriv;

		//from the left and so is positive
		//even though derivative is wrt right variable
		face = -1;
		fluxDeriv = currentDiv->getInteriorFluxDerivTerm_x(face);
		Dflux_x_right[node.interRight] += fluxDeriv;
#endif

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<init->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++init;
    }//end init loop

  //====================no flow===========================
  //mwf could try and shift over these somehow
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      if (node.isLocal)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end noit loop

}//end function

//============== 2d ======================
//this definitely won't work for multiple variables
//haven't looked at this one yet
void ConstBCECDM::applyBoundaryConditionDerivativesECDM(Petsc::StencilMM& node,
							Vec& Dflux_x_center, 
							Vec& Dflux_x_right,
							Vec& Dflux_x_back, 
							Vec& Dflux_x_front,
							Vec& Dflux_x_rightBack, 
							Vec& Dflux_x_rightFront,
							Vec& Dflux_y_center, 
							Vec& Dflux_y_back, 
							Vec& Dflux_y_right, 
							Vec& Dflux_y_left,
							Vec& Dflux_y_backRight, 
							Vec& Dflux_y_backLeft)

{
  
  //====================dirichlet===========================
  dit = dirichlet.begin();
  //std::cout<<dit<<'\t'<<dirichlet.end()<<std::endl<<std::flush;
  while(dit != dirichlet.end()) 
    {
      node.globalNode(dit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(dit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(dit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		 Dflux_x_back[node.interLeft] = 0.0; 
		 Dflux_x_front[node.interLeft] = 0.0;
		 Dflux_x_rightBack[node.interLeft] = 0.0; 
		 Dflux_x_rightFront[node.interLeft] = 0.0;
		 //
		 Dflux_y_back[node.interFront] = 0.0; 
		 Dflux_y_center[node.interFront] = 0.0;
		 Dflux_y_left[node.interFront] = 0.0; 
		 Dflux_y_right[node.interFront]= 0.0;
		 Dflux_y_backRight[node.interFront] = 0.0; 
		 Dflux_y_backLeft[node.interFront] = 0.0;
		 //
		 Dflux_y_back[node.interBack] = 0.0; 
		 Dflux_y_center[node.interBack] = 0.0;
		 Dflux_y_left[node.interBack] = 0.0; 
		 Dflux_y_right[node.interBack]= 0.0;
		 Dflux_y_backRight[node.interBack] = 0.0; 
		 Dflux_y_backLeft[node.interBack] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		//
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		
		break;
	      }
	    case FRONT:
	      {
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    case BACK:
	      {
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<dit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*dit->var)(node.center)<<'\t'<<dit->value<<endl;
      ++dit;
    } //end dit loop

  //==================== robbins ===========================
  ruit = robbins.begin();
  //std::cout<<ruit<<'\t'<<robbins.end()<<std::endl<<std::flush;
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      //mwf only set if node is local to this processor?
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(ruit->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(ruit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		 Dflux_x_back[node.interLeft] = 0.0; 
		 Dflux_x_front[node.interLeft] = 0.0;
		 Dflux_x_rightBack[node.interLeft] = 0.0; 
		 Dflux_x_rightFront[node.interLeft] = 0.0;
		 //
		 Dflux_y_back[node.interFront] = 0.0; 
		 Dflux_y_center[node.interFront] = 0.0;
		 Dflux_y_left[node.interFront] = 0.0; 
		 Dflux_y_right[node.interFront]= 0.0;
		 Dflux_y_backRight[node.interFront] = 0.0; 
		 Dflux_y_backLeft[node.interFront] = 0.0;
		 //
		 Dflux_y_back[node.interBack] = 0.0; 
		 Dflux_y_center[node.interBack] = 0.0;
		 Dflux_y_left[node.interBack] = 0.0; 
		 Dflux_y_right[node.interBack]= 0.0;
		 Dflux_y_backRight[node.interBack] = 0.0; 
		 Dflux_y_backLeft[node.interBack] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		//
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		
		break;
	      }
	    case FRONT:
	      {
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    case BACK:
	      {
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<ruit->face<<" not implemented"<<std::endl;
		exit (1);
		break;
	      }
	    }//end switch	  
	}//end isLocal
//        cout<<"D"<<'t'<<res(node.center)<<'\t'<<(*ruit->var)(node.center)<<'\t'<<ruit->value<<endl;
      ++ruit;
    } //end ruit loop


  //====================neumann===========================
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      if (node.isLocal)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(nit->n);
	  node.localNode(localNodeNumber);
	  switch(nit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		 Dflux_x_back[node.interLeft] = 0.0; 
		 Dflux_x_front[node.interLeft] = 0.0;
		 Dflux_x_rightBack[node.interLeft] = 0.0; 
		 Dflux_x_rightFront[node.interLeft] = 0.0;
		 //
		 Dflux_y_back[node.interFront] = 0.0; 
		 Dflux_y_center[node.interFront] = 0.0;
		 Dflux_y_left[node.interFront] = 0.0; 
		 Dflux_y_right[node.interFront]= 0.0;
		 Dflux_y_backRight[node.interFront] = 0.0; 
		 Dflux_y_backLeft[node.interFront] = 0.0;
		 //
		 Dflux_y_back[node.interBack] = 0.0; 
		 Dflux_y_center[node.interBack] = 0.0;
		 Dflux_y_left[node.interBack] = 0.0; 
		 Dflux_y_right[node.interBack]= 0.0;
		 Dflux_y_backRight[node.interBack] = 0.0; 
		 Dflux_y_backLeft[node.interBack] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		//
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		
		break;
	      }
	    case FRONT:
	      {
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    case BACK:
	      {
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		//
		//mwf now try and account for the fact that I switched
		//mwf flux locations?
//  		Dflux_y_back[node.interBack] = Dflux_y_back[node.interFront]; 
//  		Dflux_y_center[node.interBack] = Dflux_y_center[node.interFront];
//  		Dflux_y_left[node.interBack] = Dflux_y_left[node.interFront]; 
//  		Dflux_y_right[node.interBack]= Dflux_y_right[node.interFront];
//  		Dflux_y_backRight[node.interBack] = 
//  		  Dflux_y_backRight[node.interFront]; 
//  		Dflux_y_backLeft[node.interBack] = 
//  		  Dflux_y_backLeft[node.interFront];
		//
//  		Dflux_y_back[node.interFront] = 0.0; 
//  		Dflux_y_center[node.interFront] = 0.0;
//  		Dflux_y_left[node.interFront] = 0.0; 
//  		Dflux_y_right[node.interFront]= 0.0;
//  		Dflux_y_backRight[node.interFront] = 0.0; 
//  		Dflux_y_backLeft[node.interFront] = 0.0;
		//mwf end changes in code for new res
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<nit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++nit;
    }//end nit loop

  //====================== dummy =========================
  //zero everything
  duit = dummy.begin();
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      if (node.isLocal)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(duit->n);
	  node.localNode(localNodeNumber);
	  switch(duit->face)
	    {
	    case DUMMY_FACE:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		//
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		//
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<duit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++duit;
    } //end duit loop

  //==================== interior ====================
  //zero transverse
  //and add in extra term for interface coefficient dependence
  //on left/right front/back neighbor
  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf now have to add isLocal because I'm keeping bc's for 
      //mwf ghost region
      if (node.isLocal)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		//
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;

		//these won't work for now no div defined		
#ifdef USE_CONT_INTERFACE_TERMS_ECDM
		assert(currentDiv);
		//now get extra term
		real fluxDeriv(0.0);
		int face(0);
		//interface coefficient is upwind(Kr_{left},Kr_{right})*
		//1/2(rho_{left}+rho_{right})
	      
		//flux is from the right and so gets a negative
		//even though derivative is wrt left variable
		face = 1;
		fluxDeriv = currentDiv->getInteriorFluxDerivTerm_x(face);
		Dflux_x_center[node.interLeft] += -fluxDeriv;

		//from the left and so is positive
		//even though derivative is wrt right variable
		face = -1;
		fluxDeriv = currentDiv->getInteriorFluxDerivTerm_x(face);
		Dflux_x_right[node.interRight] += fluxDeriv;
#endif

		break;
	      }
	    case INTERIOR_Y:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;

#ifdef USE_CONT_INTERFACE_TERMS_ECDM
		assert(currentDiv);
		//now get extra term
		real fluxDeriv(0.0);
		int face(0);
		//interface coefficient is upwind(Kr_{front},Kr_{back})*
		//1/2(rho_{front}+rho_{back})
	      
		//flux is from the back and so gets a negative
		//even though derivative is wrt front variable
		face = 1;
		fluxDeriv = currentDiv->getInteriorFluxDerivTerm_y(face);
		Dflux_y_center[node.interFront] += -fluxDeriv;

		//from the front and so is positive
		//even though derivative is wrt back variable
		face = -1;
		fluxDeriv = currentDiv->getInteriorFluxDerivTerm_y(face);
		Dflux_y_back[node.interBack] += fluxDeriv;

#endif

		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<init->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++init;
    }//end init loop

  //====================no flow===========================
  //mwf could try and shift over these somehow
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      if (node.isLocal)
	{
	  //mwf now have to add isLocal because I'm keeping bc's for 
	  //mwf ghost region
	  int localNodeNumber = node.globalToLocal(noit->n);
	  node.localNode(localNodeNumber);
	  switch(noit->face)
	    {
	    case LEFT:
	      {
		 Dflux_x_center[node.interLeft] = 0.0; 
		 Dflux_x_right[node.interLeft] = 0.0;
		 Dflux_x_back[node.interLeft] = 0.0; 
		 Dflux_x_front[node.interLeft] = 0.0;
		 Dflux_x_rightBack[node.interLeft] = 0.0; 
		 Dflux_x_rightFront[node.interLeft] = 0.0;
		 //
		 Dflux_y_back[node.interFront] = 0.0; 
		 Dflux_y_center[node.interFront] = 0.0;
		 Dflux_y_left[node.interFront] = 0.0; 
		 Dflux_y_right[node.interFront]= 0.0;
		 Dflux_y_backRight[node.interFront] = 0.0; 
		 Dflux_y_backLeft[node.interFront] = 0.0;
		 //
		 Dflux_y_back[node.interBack] = 0.0; 
		 Dflux_y_center[node.interBack] = 0.0;
		 Dflux_y_left[node.interBack] = 0.0; 
		 Dflux_y_right[node.interBack]= 0.0;
		 Dflux_y_backRight[node.interBack] = 0.0; 
		 Dflux_y_backLeft[node.interBack] = 0.0;

		break;
	      }
	    case RIGHT:
	      {
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		//
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		
		break;
	      }
	    case FRONT:
	      {
		Dflux_y_back[node.interFront] = 0.0; 
		Dflux_y_center[node.interFront] = 0.0;
		Dflux_y_left[node.interFront] = 0.0; 
		Dflux_y_right[node.interFront]= 0.0;
		Dflux_y_backRight[node.interFront] = 0.0; 
		Dflux_y_backLeft[node.interFront] = 0.0;
		
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    case BACK:
	      {
		Dflux_y_back[node.interBack] = 0.0; 
		Dflux_y_center[node.interBack] = 0.0;
		Dflux_y_left[node.interBack] = 0.0; 
		Dflux_y_right[node.interBack]= 0.0;
		Dflux_y_backRight[node.interBack] = 0.0; 
		Dflux_y_backLeft[node.interBack] = 0.0;
		//
		Dflux_x_center[node.interRight] = 0.0; 
		Dflux_x_right[node.interRight] = 0.0;
		Dflux_x_back[node.interRight] = 0.0; 
		Dflux_x_front[node.interRight] = 0.0;
		Dflux_x_rightBack[node.interRight] = 0.0; 
		Dflux_x_rightFront[node.interRight] = 0.0;
		//
		Dflux_x_center[node.interLeft] = 0.0; 
		Dflux_x_right[node.interLeft] = 0.0;
		Dflux_x_back[node.interLeft] = 0.0; 
		Dflux_x_front[node.interLeft] = 0.0;
		Dflux_x_rightBack[node.interLeft] = 0.0; 
		Dflux_x_rightFront[node.interLeft] = 0.0;
		
		break;
	      }
	    default:
	      {
		std::cerr<<"ConstBCECDM.h: Face Type "<<noit->face<<" not implemented"<<std::endl;
		exit (1);
	      }
	    }//end switch
	}//end if isLocal
      ++noit;
    }//end noit loop

}//end function

//==================== interface terms
//zero derivative at interface
void  
ConstBCECDM::applyInteriorConditionToConstRelationDeriv(Petsc::StencilMM& 
							node,
							Vec& DKr,
							Vec& DRho)
{
  init = interior.begin();
  //std::cout<<init<<'\t'<<interior.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {

		DRho[node.center] = 0.0;
		DKr[node.center]  = 0.0;

		break; 
	      }
	    case INTERIOR_Y:
	      {
		DRho[node.center] = 0.0;
		DKr[node.center]  = 0.0;
		
		break;
	      }
	    default:
	      {
		cerr<<"wrong face = "<<init->face
		    <<" in applyInteriorConditionToConstRelation "<<endl;
		assert(0);
		break;
	      }

	    }//end switch
	}//end if local
      
      ++init;
    }//end while
}//end function


//try to put value into interface version of rel perms at interior boundary
//========== 1d ==========
void 
ConstBCECDM::applyInteriorConditionToInterfaceValues(Petsc::StencilMM& 
						     node,
						     const Vec* Kr,
						     const Vec* Rho,
						     Vec& KrX,Vec& RhoX)
{
  init = interior.begin();
  //std::cout<<init<<'\t'<<interior.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		//Kr and Rho should have right value loaded into center term
		KrX[node.interRightGhost]=Kr[0][node.center];
		KrX[node.interLeftGhost] =Kr[0][node.center];
		RhoX[node.interRightGhost]=Rho[0][node.center];
		RhoX[node.interLeftGhost] =Rho[0][node.center];

		break; 
	      }
	    default:
	      {
		cerr<<"wrong face = "<<init->face
		    <<" in applyInteriorConditionToConstRelation "<<endl;
		assert(0);
		break;
	      }

	    }//end switch
	}//end if local
      
      ++init;
    }//end while
}//end function

//try to put value into interface version of rel perms at interior boundary
//========== 2d ==========
void 
ConstBCECDM::applyInteriorConditionToInterfaceValues(Petsc::StencilMM& 
						     node,
						     const Vec* Kr,
						     const Vec* Rho,
						     Vec& KrX,Vec& KrY,
						     Vec& RhoX,Vec& RhoY)
{
  init = interior.begin();
  //std::cout<<init<<'\t'<<interior.end()<<std::endl<<std::flush;
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      //mwf only set if node is local to this processor?
      if (node.isInGhostRegion)
	{
	  int localNodeNumber = node.globalToLocal(init->n);
	  //mwf this needs to set flux indexes now
	  node.localNode(localNodeNumber);
	  switch(init->face)
	    {
	    case INTERIOR_X:
	      {
		//Kr and Rho should have right value loaded into center term
		KrX[node.interRightGhost]=Kr[0][node.center];
		KrX[node.interLeftGhost] =Kr[0][node.center];
		RhoX[node.interRightGhost]=Rho[0][node.center];
		RhoX[node.interLeftGhost] =Rho[0][node.center];

		break; 
	      }
	    case INTERIOR_Y:
	      {
		//Kr and Rho should have right value loaded into center term
		KrY[node.interBackGhost]  =Kr[0][node.center];
		KrY[node.interFrontGhost] =Kr[0][node.center];
		RhoY[node.interBackGhost] =Rho[0][node.center];
		RhoY[node.interFrontGhost]=Rho[0][node.center];
	
		break;
	      }
	    default:
	      {
		cerr<<"wrong face = "<<init->face
		    <<" in applyInteriorConditionToConstRelation "<<endl;
		assert(0);
		break;
	      }

	    }//end switch
	}//end if local
      
      ++init;
    }//end while
}//end function

////////////////////////////
void  ConstBCECDM::adjustTimeIntegrationTolerance(Petsc::StencilMM& node, 
				    Vec& tol, real newTol)
{
  //set larger tolerances on neumann,robbins interior and dummy boundaries
  nit = neumann.begin();
  while(nit != neumann.end()) 
    {
      node.globalNode(nit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  tol(globalCenter) = newTol;
	}
      ++nit;
    }
  ruit = robbins.begin();
  while(ruit != robbins.end()) 
    {
      node.globalNode(ruit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  tol(globalCenter) = newTol;
	}
      ++ruit;
    }
  //
  duit = dummy.begin();
  while(duit != dummy.end()) 
    {
      node.globalNode(duit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  tol(globalCenter) = newTol;
	}
      ++duit;
    }
  //
  init = interior.begin();
  while(init != interior.end()) 
    {
      node.globalNode(init->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  tol(globalCenter) = newTol;
	}
      ++init;
    }
  noit = noflow.begin();
  while(noit != noflow.end()) 
    {
      node.globalNode(noit->n);
      int globalCenter = node.center;
      if (node.isLocal)
	{
	  tol(globalCenter) = newTol;
	}
      ++noit;
    }

}

}//Daetk




