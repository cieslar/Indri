#include<DMContainer.hpp>
#include<chealpix.h>
#include<fstream>
#include<cstring>
#include<NE2001.hpp>
#include<cmath>
#include<NE2001.hpp>
#include<ConfigContainer.hpp>

using namespace std;

RealType GalLDeg2HEALPixPhi(const RealType fLDeg);
RealType HEALPixPhi2GalLDeg(const RealType fPhi);
RealType GalBDeg2HEALPixTheta(const RealType fBDeg);
RealType HEALPixTheta2GalBDeg(const RealType fTheta);

RealType GalLDeg2HEALPixPhi(const RealType fLDeg)
{
	return fLDeg*M_PI/180.;
}

RealType HEALPixPhi2GalLDeg(const RealType fPhi)
{
	return fPhi*180./M_PI;
}

RealType GalBDeg2HEALPixTheta(const RealType fBDeg)
{
	return (90.-fBDeg)*M_PI/180.;
}

RealType HEALPixTheta2GalBDeg(const RealType fTheta)
{
	return   (M_PI/2.-fTheta)*180./M_PI;
}

/////////////////////////////////////////////

CDMContainer::CDMContainer(const char * FileName)
{
	ReadBinaryGeometry(FileName);
	m_PrepareContainer(m_nSide, m_nDistNum, m_fMaxDist);
	ReadBinary(FileName);
}

CDMContainer::CDMContainer(const long int nSide, const long int nDistNum, const RealType fMaxDist)
{
	m_PrepareContainer(nSide, nDistNum, fMaxDist);
}

void CDMContainer::m_PrepareContainer(const long int nSide, const long int nDistNum, const RealType fMaxDist)
{
	
	m_nSide= nSide;
	m_nPix  = nside2npix(m_nSide);
	m_nDistNum = nDistNum;
	m_fMaxDist = fMaxDist;
	m_fDDist = m_fMaxDist/m_nDistNum;
	

	m_ppDM = new float*[m_nDistNum];

	for(int d(0);d<m_nDistNum;++d)
	{
		m_ppDM[d]=new float[m_nPix];
	}

	unsigned long long int nZlicz(0);
	for(int d(0);d<m_nDistNum;++d) for(int i(0);i<m_nPix;++i)
	{
		m_ppDM[d][i]=0.;
		nZlicz++;
	}

}




void CDMContainer::WriteSliceBinary(const long int nSlice, const char* FileName)
{
	size_t nSliceSize ( m_nPix*sizeof(float) );
	if(nSlice<0 || nSlice>m_nDistNum) ERROR("Outside of slices' range");

	std::ofstream Out(FileName,std::ofstream::binary);

	if (!Out.is_open()) ERROR("Can't open file for writing");

	Out.write(reinterpret_cast<char *>(m_ppDM[nSlice]),nSliceSize);

	Out.close();
}

void CDMContainer::ReadSliceBinary(const long int nSlice, const char* FileName)
{
	size_t nSliceSize ( m_nPix*sizeof(float) );
	if(nSlice<0 || nSlice>m_nDistNum) ERROR("Outside of slices' range");

	std::ifstream In(FileName,std::ifstream::binary);

	
	In.seekg( 0, std::ios::end );
	size_t nFileSize =  In.tellg();
	In.seekg (0, In.beg);

	if(nFileSize!=nSliceSize) ERROR("DMSliceContainer size mismatch. Rading not possible.");

	In.read(reinterpret_cast<char *>(m_ppDM[nSlice]),nSliceSize);
	In.close();
}





void CDMContainer::WriteBinary(const char* FileName)
{

	std::ofstream Out(FileName,std::ofstream::binary);

	if (!Out.is_open()) ERROR("Can't open file for writing");

	Out.write(reinterpret_cast<char *>(&m_nSide), sizeof(m_nSide) );
	Out.write(reinterpret_cast<char *>(&m_nDistNum), sizeof(m_nDistNum) );
	Out.write(reinterpret_cast<char *>(&m_fMaxDist), sizeof(m_fMaxDist) );



	size_t nSliceSize ( m_nPix*sizeof(float) );
	for(int i(0);i<m_nDistNum;++i)
	{
		Out.write(reinterpret_cast<char *>(m_ppDM[i]),nSliceSize);
	}

	Out.close();
}


void CDMContainer::ReadBinaryGeometry(const char* FileName)
{
	std::ifstream In(FileName,std::ifstream::binary);

	In.read(reinterpret_cast<char *>(&m_nSide), sizeof(m_nSide) );
	In.read(reinterpret_cast<char *>(&m_nDistNum), sizeof(m_nDistNum) );
	In.read(reinterpret_cast<char *>(&m_fMaxDist), sizeof(m_fMaxDist) );
	
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();

	pCfgPtr->Sim.nHEALPixNSide=m_nSide;
	pCfgPtr->Sim.nDMDistSlices=m_nDistNum;
	pCfgPtr->Sim.fMaxDist=m_fMaxDist;
	m_nPix  = nside2npix(m_nSide);
	m_fDDist = m_fMaxDist/m_nDistNum;

	In.close();
}



void CDMContainer::ReadBinary(const char* FileName)
{
	std::ifstream In(FileName,std::ifstream::binary);

	size_t nSliceSize ( m_nPix*sizeof(float) );
	size_t nSize (sizeof(m_nSide)+sizeof(m_nDistNum)+sizeof(m_fMaxDist)
                     +m_nDistNum*nSliceSize);
	
	In.seekg( 0, std::ios::end );
	size_t nFileSize =  In.tellg();
	In.seekg (0, In.beg);

	if(nFileSize!=nSize) ERROR("DMContainer size mismatch. Rading not possible.");

	In.seekg(sizeof(m_nSide)+sizeof(m_nDistNum)+sizeof(m_fMaxDist));

	for(int i(0);i<m_nDistNum;++i)
	{
		In.read(reinterpret_cast<char *>(m_ppDM[i]),nSliceSize);
	}
	In.close();
}






CDMContainer::~CDMContainer()
{
	for(int i(0);i<m_nDistNum;++i)
	{
		delete m_ppDM[i];
	}
	delete m_ppDM;
}

//TODO dodać interpolację
RealType CDMContainer::GetDM(const RealType fLDeg, const RealType fBDeg, const RealType fDist) const
{
	long int nPix;
	long int nDist(static_cast<long int>(fDist/m_fDDist)-1 );

	if (nDist<0) nDist=0;
	if(nDist>m_nDistNum) 
	{
		cerr<<nDist<<" "<<m_nDistNum<<" "<<fLDeg<<" "<<fBDeg<<" "<<fDist<<endl;
		WARNING("Requested DM oustide of computed space");
	}

	ang2pix_ring(m_nSide,GalBDeg2HEALPixTheta(fBDeg), GalLDeg2HEALPixPhi(fLDeg), &nPix);
	return m_ppDM[nDist][nPix];

}


void CDMContainer::ComputeTheData()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	RealType fTheta(0.), fPhi(0.), fDist(0.);
	for(long int nDist(0);nDist<m_nDistNum;++nDist) 
	{
		if(pCfgPtr->Sim.bVerbose) cout<<"Do konca pozostało map: "<<m_nDistNum-nDist<<endl;
		ComputeSlice(nDist);
	}
}


void CDMContainer::ComputeSlice(const long int nSlice)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	RealType fProc(0.),fAvgTime(0.);
	int nWypisz(1);
	TTimestamp Beg(0), DT(0), WorstDT(0), BestDT(0);
	RealType fTheta(0.), fPhi(0.), fDist(0.);

	fDist = (nSlice+1)*m_fDDist;
	if(pCfgPtr->Sim.bVerbose ) cout<<m_nPix<<" pixels to compute at distance: "<<fDist<<"[kpc]"<<endl;

	for(long int nPix(0);nPix<m_nPix;++nPix)
	{
		if(pCfgPtr->Sim.bVerbose ) 
		{
			fProc=100.*static_cast<float>(nPix+1)/m_nPix;
			if( fProc>=nWypisz) 
			{
				nWypisz+=1;
				cout<<"Slice done in "<<static_cast<int>(fProc)<<"%. ";
				cout<<"Computation times till now: Worst: "<<static_cast<float>(WorstDT)/1000000.<<"[s] "
							         <<"Best: "<<static_cast<float>(BestDT)/1000000. <<"[s] "
							      <<"Average: "<<fAvgTime/1000000.<<"[s] "<<endl;
			}
		}
		Beg=GetTimestamp();

		pix2ang_ring(m_nSide,nPix,&fTheta, &fPhi);
		m_ppDM[nSlice][nPix] = fort_dmdsm(HEALPixPhi2GalLDeg(fPhi)*M_PI/180., HEALPixTheta2GalBDeg(fTheta)*M_PI/180., -1, 0., fDist);
//RADIANY!!! TODO
		DT=GetTimestamp() - Beg;
		if(DT>WorstDT) WorstDT=DT;
		if(DT<BestDT) BestDT=DT;
		if(nPix==0) BestDT=DT;
		fAvgTime= fAvgTime* nPix/(nPix+1) + static_cast<float>(DT)/(nPix+1);
			
	}
}


















#if 0
//funkcja testowa

int main()
{
	long int nNSide(128);
	long int nNPix;//number of piksel
	nNPix=nside2npix(nNSide);

	cout<<"nNSide: "<<nNSide<<" nNPix: "<<nNPix<<endl;
	long int Pix;
	RealType fL(0.),fB(-90.);
	cout<<fB<<" "<<GalBDeg2HEALPixTheta(fB)*180./M_PI<<" "<<HEALPixTheta2GalBDeg(GalBDeg2HEALPixTheta(fB))<<endl;
	cout<<fL<<" "<<GalLDeg2HEALPixPhi(fL)*180./M_PI<<" "<<HEALPixPhi2GalLDeg(GalLDeg2HEALPixPhi(fL))<<endl;
	
	ang2pix_ring(nNSide,GalBDeg2HEALPixTheta(fB), GalLDeg2HEALPixPhi(fL), &Pix);
	cout<<nNSide<<" "<<fB<<" "<<fL<<" "<<Pix<<endl;

	RealType fTheta, fPhi;


	long npix, nside ; 
	char order[10] ; 
	npix= get_fits_size("pixel_coords_map_ring_galactic_res7.fits", &nside, order); 
	cout<<npix<<endl;
	for(int i(0);i<nNPix;++i)
	{
		pix2ang_ring(nNSide,i,&fTheta, &fPhi);
		cout<<nNSide<<" "<<i<<" "<<HEALPixPhi2GalLDeg(fPhi)<<" "<<HEALPixTheta2GalBDeg(fTheta)<<endl;
	}
	long int nside;
	char coord,ordering;
	float *pData=read_healpix_map("pixel_coords_map_ring_galactic_res7.fits",&nside,&coord,&ordering);
	cout<<nside<<" "<<coord<<" "<<ordering<<endl;
npix=nside2npix(nside);
	for(int i(0);i<npix;++i)
	{
		pix2ang_ring(nside,i,&fTheta, &fPhi);
		fPhi=HEALPixPhi2GalLDeg(fPhi);
		cout<<i<<" "<<pData[i]-fPhi<<" "<<pData[i]<<" "<<fPhi<<endl;
	}
	return 1;
}

#endif
