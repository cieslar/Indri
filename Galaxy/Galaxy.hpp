#pragma once

#include<Definitions.hpp>
#include<Coordinates.hpp>
#include<Randomizer.hpp>





using namespace std;

///An enumerator to keep the galactic arm's names.
enum GalacticArmEnu {NormaArm=1, CarinaSagittariusArm=2, PerseusArm=3, CruxScutumArm=4};

const RealType gfG=4.49824e-2; ///< [kpc^3 (1e10Msun)^-1 Myr^-2]

///Class for the computation of the galactic gravitational potential, acceleration of the pulsar and initial position within galactic arms.
class CGalaxy
{
protected:
	CRandomizer* m_pRanPtr;

	RealType  m_SunOmega;///< Angular velocity of the Sun [rad]
	//const Cartesian m_SunCurrentPosition(0., 8.5, 0.);///< Current sun's position int the cartesian coordinates in [kpc]
	Cartesian m_SunCurrentPosition;///< Current sun's position int the cartesian coordinates in [kpc]

	///  Computes the m_SunOmega
	RealType m_ComputeSunOmega() ;///< [rad]

	///  Computes the angle that the Galaxy has rotated back to the fTimePoint assuming rigid rotation
	RealType m_ComputeRigidRotationAngle(const RealType fTimePoint) const;///< [rad]


	///The potential's components [kpc^2 Myr^-2]
	RealType m_MiyamotoNagaiDisk(const RealType fMass, const RealType fParA, const RealType fParB, const Cartesian Position) const;
	RealType m_MiyamotoNagaiBuldge(const RealType fMass, const RealType fParB, const Cartesian Position) const;
	RealType m_Hernquist(const RealType fMass, const RealType fParC, const Cartesian Position) const;
	RealType m_Paczynski(const RealType fMass, const RealType fRCore, const RealType fRCut, const Cartesian Position) const;
	
	///The acceleration's components [kpc Myr^-2]
	Cartesian m_DivMiyamotoNagaiDisk(const RealType fMass, const RealType fParA, const RealType fParB, const Cartesian Position) const;
	Cartesian m_DivMiyamotoNagaiBuldge(const RealType fMass, const RealType fParB, const Cartesian Position) const;
	Cartesian m_DivHernquist(const RealType fMass, const RealType fParC, const Cartesian Position) const;
	Cartesian m_DivPaczynski(const RealType fMass, const RealType fRCore, const RealType fRCut, const Cartesian Position) const;


	///Initial positions and velocity functions for the pulsars
	RealType m_YusifovKucukRhoDist(const RealType fRho) const;
	RealType m_RandomPulsarRho() const;
	RealType m_SpiralArmTheta(const RealType fRho, const int nArmNum) const;
	RealType m_PaczysnkiKicksDist(const RealType fVelocity) const;


	RealType fR;
	RealType fRho;

	void m_ComputePartialElements(const Cartesian Position);

public:

	CGalaxy();
	~CGalaxy();

	RealType GalacticPotential(const Cartesian Position) ;///< [kpc^2/Myr^2]
	Cartesian GalacticAcceleration(const Cartesian Position) ;///< [kpc/Myr^2] for each coordinate
	Cartesian RandomPulsarInitPos() const;///< [kpc] for each coordinate
	Cartesian PaczynskiKick() const;///< [kpc/Myr] for each coordinate
	Cartesian HobbsKick() const;///<[kpc/Myr]

};


class CTestGalaxy : public CGalaxy
{
public:
	CTestGalaxy();
	~CTestGalaxy(){};

	void ArmsTest();
	void KickVelocityTest();
	void PositionTest();
	void CoordinatesTransformTest();
	void PotentialTestDirection();
	void PotentialTest();

	void PaczynksiKickTest(const int nNum);
	void HobbsKickTest(const int nNum);

	void RunAll();
};


