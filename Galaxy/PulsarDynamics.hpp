#pragma once

#include<Galaxy.hpp>
//#include<ConfigContainer.hpp>
#include<Definitions.hpp>
#include<BasicStar.hpp>
#include<Coordinates.hpp>

class CPulsarDynamics: public virtual CBasicStar, public CGalaxy
{
protected:
	void m_EvolveRK3Adaptive(const RealType fMaxTimeStep);//w [Myr]
	void m_EvolveNewtonForward(const RealType fTimeStep);//w [Myr]
	void m_EvolveLeapfrog(const RealType fTimeStep);//w [Myr]
	void m_VelocityHalfLeap(const RealType fTimeStep);
	void m_EvolvePositionVerlet(const RealType fTimeStep);
	void m_EvolvePositionVerletDynamic(const RealType fMaxTimeStep);

	Galactic m_Project(const RealType fProjectionTime=0.);//w Myr czas kiedy robie projekcje 0.=teraz

public:

	CPulsarDynamics(){Velocity.fX=Velocity.fY=Velocity.fZ=0.;};
	~CPulsarDynamics(){};
	Cartesian GetPosition() const;
	void Init(const RealType fPositionStartTime=0.);
	void Evolve(const RealType fMaxTimeStep);
	void EvolveEnd(const RealType fMaxTimeStep);
	void EvolveInit(const RealType fMaxTimeStep);

	//void ProjectTest();
};

