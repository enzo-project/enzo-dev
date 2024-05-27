/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/  modified2:  John Wise, Greg Bryan, Britton Smith, Cameron Hummels,
/              Matt Turk
/  date:       May, 2011 (converting from Star to ActiveParticle)
/
/  PURPOSE:
/
************************************************************************/

#ifndef __PARTICLE_ATTRIBUTE_HANDLER_H
#define __PARTICLE_ATTRIBUTE_HANDLER_H

/* http://www.gamedev.net/topic/474803-c-template-pointer-to-member/ */
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <hdf5.h>
class ActiveParticleType;

class ParticleAttributeHandler
{

  public:

    std::string name;
#ifdef USE_MPI
    MPI_Datatype mpitype;
#endif
    hid_t hdf5type;
    int element_size;
    int offset;

    virtual void SetAttribute(char **buffer, ActiveParticleType *pp) = 0;

    virtual int GetAttribute(char **buffer, ActiveParticleType *pp)  = 0;

    virtual void PrintAttribute(ActiveParticleType *pp)  = 0;

};

template <class APClass, typename Type, Type APClass::*var>
class Handler : public ParticleAttributeHandler
{
  public:

    Handler(std::string name, int offset = 0) {
        this->name = name;
        this->offset = offset;

        /* Can't use a switch */
        if (typeid(Type) == typeid(int)) {
#ifdef USE_MPI
            this->mpitype = IntDataType;
#endif
            this->hdf5type = HDF5_INT;
        } else if (typeid(Type) == typeid(float)) {
#ifdef USE_MPI
            this->mpitype = FloatDataType;
#endif
            this->hdf5type = HDF5_REAL;
        } else if (typeid(Type) == typeid(double)) {
#ifdef USE_MPI
            this->mpitype = MPI_DOUBLE;
#endif
            this->hdf5type = HDF5_R8;
        } else if (typeid(Type) == typeid(FLOAT)) {
#ifdef USE_MPI
            this->mpitype = FLOATDataType;
#endif
            this->hdf5type = HDF5_PREC;
        } else {
            ENZO_FAIL("Unrecognized data type");
        }
        this->element_size = sizeof(Type);
    }

    void SetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
        pp->*var = *(pb++);
        *buffer = (char *) pb;
    }

    int GetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
        *(pb++) = pp->*var;
        *buffer = (char *) pb;
        return this->element_size;
    }

    void PrintAttribute(ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        std::cout << std::setprecision(15) << this->name << ": " << pp->*var;
    }

};

template <class APClass, typename Type, int N, Type (APClass::*var)[N]>
class ArrayHandler : public ParticleAttributeHandler
{
  public:

    ArrayHandler(std::string name, int offset = 0) {
        this->name = name;
        this->offset = offset;

        /* Can't use a switch */
        if (typeid(Type) == typeid(int)) {
#ifdef USE_MPI
            this->mpitype = IntDataType;
#endif
            this->hdf5type = HDF5_INT;
        } else if (typeid(Type) == typeid(float)) {
#ifdef USE_MPI
            this->mpitype = FloatDataType;
#endif
            this->hdf5type = HDF5_REAL;
        } else if (typeid(Type) == typeid(double)) {
#ifdef USE_MPI
            this->mpitype = MPI_DOUBLE;
#endif
            this->hdf5type = HDF5_R8;
        } else if (typeid(Type) == typeid(FLOAT)) {
#ifdef USE_MPI
            this->mpitype = FLOATDataType;
#endif
            this->hdf5type = HDF5_PREC;
        } else {
            ENZO_FAIL("Unrecognized data type");
        }
	const char *_name = this->name.c_str();
	/* particle_position and particle_velocity are not actually stored as arrays */
	if(strncmp(_name, "particle_", 9) == 0)
	  this->element_size = sizeof(Type);  
	else
	  this->element_size = sizeof(Type)*N;
    }

    void SetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);
	
	/* For everything except particle_position and particle_velocity 
	 * we need to loop over the size of the arrays
	 */
	const char *_name = this->name.c_str();
	
	if(strncmp(_name, "particle_", 9) == 0) {
	  (pp->*var)[this->offset] = *(pb++);
	  *buffer = (char *) pb;
	}
	else {
	  for(int ii = 0; ii < N; ii++) {
	    (pp->*var)[ii] =  *(pb++);
	    *(buffer)  = (char *) pb;
	  }
	}
    }

    int GetAttribute(char **buffer, ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        Type *pb = (Type *)(*buffer);

	/* For everything except particle_position and particle_velocity 
	 * we need to loop over the size of the arrays
	 */
	const char *_name = this->name.c_str();
	if(strncmp(_name, "particle_", 9) == 0) {
	  *(pb++) = (pp->*var)[this->offset];
	  *buffer = (char *) pb;
	}
	else {
	  for(int ii = 0; ii < N; ii++) {
	    *(pb++) = (pp->*var)[ii];
	    *(buffer) = (char *) pb;
	  }
	}
        return this->element_size;
    }

    void PrintAttribute(ActiveParticleType *pp_) {
        APClass *pp = static_cast<APClass*>(pp_);
        std::cout << this->name << ": " << (pp->*var)[this->offset];
    }

};

typedef std::vector<ParticleAttributeHandler*> AttributeVector ;

#endif
