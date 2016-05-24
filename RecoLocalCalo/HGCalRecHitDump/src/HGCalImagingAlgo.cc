#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalImagingAlgo.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"

// Return a vector of clusters from a collection of EcalRecHits:
std::vector<reco::BasicCluster> HGCalImagingAlgo::makeClusters(
                                  const HGCRecHitCollection& hits)
{

  const HGCalDDDConstants* ddd = &(geometry->topology().dddConstants());


  current_v.clear();
  clusters_v.clear();
  cluster_offset = 0;

  std::vector<std::vector<Hexel> > points(2*(maxlayer+1)); //a vector of vectors of hexels, one for each layer
  //@@EM todo: the number of layers should be obtained programmatically - the range is 1-n instead of 0-n-1...


  if (verbosity < pINFO)
    {
      std::cout << "-------------------------------------------------------------" << std::endl;
      std::cout << "HGC Imaging algorithm invoked for HGEE" << std::endl;
      std::cout << "delta_c " << delta_c << " kappa " << kappa <<std::endl;
    }

  //loop over all hits and create the Hexel structure, skip energies below ecut
  for (unsigned int i=0;i<hits.size();++i) {
    const HGCRecHit& hgrh = hits[i];
    if(hgrh.energy() < ecut) continue; 
    DetId detid = hgrh.detid();

    int layer = HGCEEDetId(detid).layer()+int(HGCEEDetId(detid).zside()>0)*(maxlayer+1);

    // determine whether this is a half-hexagon
    // (copied from Lindsey's code not (yet?) available in release - is this even right ?

    bool isHalf = false;
    if(ddd!=0){
      const HGCalDetId hid(detid);
      const int waferType = ddd->waferTypeT(hid.waferType());  
      isHalf = ddd->isHalfCell(waferType,hid.cell());
    }
    points[layer].push_back(Hexel(hgrh,detid,isHalf,geometry));
  }
  //assign all hits in each layer to a cluster core or halo
  for (unsigned int i = 0; i <= 2*maxlayer+1; ++i) {
    double maxdensity = calculateLocalDensity(points[i]);
    // std::cout << "layer " << i << " max density " << maxdensity 
    // 	      << " total hits " << points[i].size() << std::endl;
    calculateDistanceToHigher(points[i]);
    findAndAssignClusters(points[i],maxdensity);
    //    std::cout << "found " << nclusters << " clusters" << std::endl;
  }
  //make the cluster vector
  reco::CaloID caloID = reco::CaloID::DET_HGCAL_ENDCAP;
  for (unsigned int i = 0; i < current_v.size(); i++){
    double energy = 0;
    Point position;
    std::vector< std::pair<DetId, float> > thisCluster;
    
    if( doSharing ) {
      std::vector<unsigned> seeds = findLocalMaximaInCluster(current_v[i]);
      std::vector<std::vector<double> > fractions;
      shareEnergy(current_v[i],seeds,fractions);
      
      for( unsigned isub = 0; isub < fractions.size(); ++isub ) {
	double effective_hits = 0.0;
	double energy  = calculateEnergyWithFraction(current_v[i],fractions[isub]);
	Point position = calculatePositionWithFraction(current_v[i],fractions[isub]);
	
	for( unsigned ihit = 0; ihit < fractions[isub].size(); ++ihit ) {
	  const double fraction = fractions[isub][ihit];
	  effective_hits += fraction;
	  thisCluster.emplace_back(current_v[i][ihit].detid,fraction);
	}
	
	if (verbosity < pINFO)
	  { 
	    std::cout << "******** NEW CLUSTER ********" << std::endl;
	    std::cout << "No. of crystals = " << effective_hits << std::endl;
	    std::cout << "     Energy     = " << energy << std::endl;
	    std::cout << "     Phi        = " << position.phi() << std::endl;
	    std::cout << "     Eta        = " << position.eta() << std::endl;
	    std::cout << "*****************************" << std::endl;
	  }
	clusters_v.push_back(reco::BasicCluster(energy, position, caloID, thisCluster, 
						reco::CaloCluster::hgcal_em));
	thisCluster.clear();
      }
    } else {
      position = calculatePosition(current_v[i]);    
      std::vector< Hexel >::iterator it;
      for (it = current_v[i].begin(); it != current_v[i].end(); it++)
	{
	  energy += (*it).weight;
	  thisCluster.emplace_back(std::pair<DetId, float>((*it).detid,((*it).isHalo?0.:1.)));
	}
      if (verbosity < pINFO)
	{ 
	  std::cout << "******** NEW CLUSTER ********" << std::endl;
	  std::cout << "No. of crystals = " << current_v.size() << std::endl;
	  std::cout << "     Energy     = " << energy << std::endl;
	  std::cout << "     Phi        = " << position.phi() << std::endl;
	  std::cout << "     Eta        = " << position.eta() << std::endl;
	  std::cout << "*****************************" << std::endl;
	}
      clusters_v.push_back(reco::BasicCluster(energy, position, caloID, thisCluster, 
					      reco::CaloCluster::hgcal_em));
    }    
  }
  return clusters_v; 
}  

math::XYZPoint HGCalImagingAlgo::calculatePosition(std::vector<Hexel> &v){
  float total_weight = 0.;
  float x = 0.;
  float y = 0.;
  float z = 0.;
  for (unsigned int i = 0; i < v.size(); i++){
    total_weight += v[i].weight;
    x += v[i].x*v[i].weight;
    y += v[i].y*v[i].weight;
    z += v[i].z*v[i].weight;
  }
  
  return math::XYZPoint( x/total_weight, 
			 y/total_weight, 
			 z/total_weight );
} 

double HGCalImagingAlgo::distance(const Hexel &pt1, const Hexel &pt2){
  const GlobalPoint position1( std::move( geometry->getPosition( pt1.detid ) ) );
  const GlobalPoint position2( std::move( geometry->getPosition( pt2.detid ) ) );
  return sqrt(pow(pt1.x - pt2.x, 2) + pow(pt1.y - pt2.y, 2));
}


double HGCalImagingAlgo::calculateLocalDensity(std::vector<Hexel> &lp){
  double maxdensity = 0.;
  for(unsigned int i = 0; i < lp.size(); i++){
    for(unsigned int j = 0; j < lp.size(); j++){
      if(distance(lp[i],lp[j]) < delta_c){
	lp[i].rho += lp[j].weight;
	if(lp[i].rho > maxdensity) maxdensity = lp[i].rho;
      }
    }
  }
  return maxdensity;
}

double HGCalImagingAlgo::calculateDistanceToHigher(std::vector<Hexel> &lp){
  

  //sort vector of Hexels by decreasing local density
  std::vector<size_t> rs = sorted_indices(lp);

  double maxdensity = 0.0;
  int nearestHigher = -1;


  if(rs.size()>0) 
    maxdensity = lp[rs[0]].rho;
  else
    return maxdensity; // there are no hits
  double dist = 50.0;
  //start by setting delta for the highest density hit to 
  //the most distant hit - this is a convention

  for(unsigned int j = 0; j < lp.size(); j++){
    double tmp = distance(lp[rs[0]], lp[j]);
    dist = tmp > dist ? tmp : dist;
  }
  lp[rs[0]].delta = dist;
  lp[rs[0]].nearestHigher = nearestHigher;

  //now we save the largest distance as a starting point
  
  double max_dist = dist;
  
  for(unsigned int oi = 1; oi < lp.size(); oi++){ // start from second-highest density
    dist = max_dist;
    unsigned int i = rs[oi];
    // we only need to check up to oi since hits 
    // are ordered by decreasing density
    // and all points coming BEFORE oi are guaranteed to have higher rho 
    // and the ones AFTER to have lower rho
    for(unsigned int oj = 0; oj < oi; oj++){ 
      unsigned int j = rs[oj];
      double tmp = distance(lp[i], lp[j]);
      if(tmp <= dist){ //this "<=" instead of "<" addresses the (rare) case when there are only two hits
	dist = tmp;
	nearestHigher = j;
      }
    }
    lp[i].delta = dist;
    lp[i].nearestHigher = nearestHigher; //this uses the original unsorted hitlist 
  }
  return maxdensity;
}

int HGCalImagingAlgo::findAndAssignClusters(std::vector<Hexel> &lp, double maxdensity){

  //this is called once per layer...
  //so when filling the cluster temporary vector of Hexels we resize each time by the number 
  //of clusters found. This is always equal to the number of cluster centers...

  unsigned int clusterIndex = 0;

  std::vector<size_t> rs = sorted_indices(lp); // indices sorted by decreasing rho
  std::vector<size_t> ds = sort_by_delta(lp); // sort in decreasing distance to higher


  for(unsigned int i =0; i < lp.size(); i++){

    //    std::cout << " delta " << lp[ds[i]].delta << " rho " << lp[ds[i]].rho << std::endl;
    if(lp[ds[i]].delta < delta_c) break; // no more cluster centers to be looked at 
    if(lp[ds[i]].rho < maxdensity/kappa || lp[ds[i]].rho<0.001) continue; 
    //skip this as a potential cluster center because it fails the density cut

    lp[ds[i]].clusterIndex = clusterIndex;
    clusterIndex++;
  }

  //at this point clusterIndex is equal to the number of cluster centers - if it is zero we are 
  //done
  if(clusterIndex==0) return clusterIndex;

  //assign to clusters, using the nearestHigher set from previous step (always set except 
  // for top density hit that is skipped...
  for(unsigned int oi =1; oi < lp.size(); oi++){
    unsigned int i = rs[oi];
    int ci = lp[i].clusterIndex;
    if(ci == -1){
      lp[i].clusterIndex =  lp[lp[i].nearestHigher].clusterIndex;
    }
  }

  //make room in the temporary cluster vector for the additional clusterIndex clusters 
  // from this layer
  current_v.resize(cluster_offset+clusterIndex);

  //assign points closer than dc to other clusters to border region
  //and find critical border density
  std::vector<double> rho_b(clusterIndex,0.);

  for(unsigned int i = 0; i < lp.size(); i++){
    int ci = lp[i].clusterIndex;
    bool flag_isolated = true;
    if(ci != -1){
      for(unsigned int j = 1; j < lp.size(); j++){
	//check if the hit is not within d_c of another cluster
	if(lp[j].clusterIndex!=-1){
	  float dist = distance(lp[j],lp[i]);
	  if(dist < delta_c && lp[j].clusterIndex!=ci){
	    //in which case we assign it to the border
	    lp[i].isBorder = true;
	    break;
	  }
	  if(dist < delta_c && lp[j].clusterIndex==ci){
	    //in which case we assign it to the border
	    flag_isolated = false;
	  }
	}
      }
      if(flag_isolated) lp[i].isBorder = true; //the hit is more than delta_c from any of its brethren
    }	  
    if(lp[i].isBorder && rho_b[ci] < lp[i].rho)
      rho_b[ci] = lp[i].rho;
  }

  //flag points in cluster with density < rho_b as halo points, then fill the cluster vector 
  for(unsigned int i = 1; i < lp.size(); i++){
    int ci = lp[i].clusterIndex;
    if(ci!=-1 && lp[i].rho < rho_b[ci])
      lp[i].isHalo = true;
    if(lp[i].clusterIndex!=-1) 
      current_v[ci+cluster_offset].push_back(lp[i]);
  }

  //prepare the offset for the next layer if there is one
  cluster_offset += clusterIndex;
  return clusterIndex;
}

// find local maxima within delta_c, marking the indices in the cluster
std::vector<unsigned>&& HGCalImagingAlgo::findLocalMaximaInCluster(const std::vector<Hexel>& cluster) {
  std::vector<unsigned> result;
  std::vector<bool> seed(cluster.size(),true);
 
  for( unsigned i = 0; i < cluster.size(); ++i ) {    
    for( unsigned j = 0; j < cluster.size(); ++j ) {
      if( distance(cluster[i],cluster[j]) < delta_c && i != j) {
	if( cluster[i].weight < cluster[j].weight ) seed[i] = false;
      }
    }
  }

  for( unsigned i = 0 ; i < cluster.size(); ++i ) {
    if( seed[i] ) result.push_back(i);
  }

  std::cout << "Found " << result.size() << " sub-clusters!" << std::endl;

  return std::move(result);
}

math::XYZPoint&& HGCalImagingAlgo::calculatePositionWithFraction(const std::vector<Hexel>& hits,
								 const std::vector<double>& fractions) {  
  double norm(0.0), x(0.0), y(0.0), z(0.0);
  for( unsigned i = 0; i < hits.size(); ++i ) {
    const double weight = fractions[i]*hits[i].weight;
    norm += weight;
    x += weight*hits[i].x;
    y += weight*hits[i].y;
    z += weight*hits[i].z;
  }
  math::XYZPoint result(x,y,z);
  double norm_inv = 1.0/norm;
  result *= norm_inv;
  return std::move(result);
}

double HGCalImagingAlgo::calculateEnergyWithFraction(const std::vector<Hexel>& hits,
						     const std::vector<double>& fractions) {
  double result = 0.0;
  for( unsigned i = 0 ; i < hits.size(); ++i ) {
    result += fractions[i]*hits[i].weight;
  }
  return result;
}

void HGCalImagingAlgo::shareEnergy(const std::vector<Hexel>& incluster,
				   const std::vector<unsigned>& seeds, 
				   std::vector<std::vector<double> >& outclusters) {
  std::vector<bool> isaseed(incluster.size(),false);
  outclusters.clear();
  outclusters.resize(seeds.size());
  std::vector<Point> centroids(seeds.size());
  std::vector<double> energies(seeds.size());

  // create quick seed lookup
  for( unsigned i = 0; i < seeds.size(); ++i ) {
    isaseed[seeds[i]] = true;
  }

  // initialize clusters to be shared
  // centroids start off at seed positions
  // seeds always have fraction 1.0, to stabilize fit
  for( unsigned i = 0; i < seeds.size(); ++i ) {
    outclusters[i].resize(incluster.size(),0.0);
    for( unsigned j = 0; j < incluster.size(); ) {
      if( j == seeds[i] ) {
	outclusters[i][j] = 1.0;
	centroids[i] = math::XYZPoint(incluster[j].x,incluster[j].y,incluster[j].z);
	energies[i]  = incluster[j].weight; 
      } 
    }
  }

  // run the fit while we are less than max iterations, and clusters are still moving
  const double minFracTot = 1e-20;
  unsigned iter = 0;
  const unsigned iterMax = 50;
  double diff = std::numeric_limits<double>::max();
  const double stoppingTolerance = 1e-8;
  const double toleranceScaling = std::pow(std::max(1.0,seeds.size()-1.0),2.0);
  std::vector<Point> prevCentroids;
  while( iter++ < iterMax && diff > stoppingTolerance*toleranceScaling ) {
    std::vector<double> frac(seeds.size()), dist2(seeds.size());
    for( unsigned i = 0; i < incluster.size(); ++i ) {
      const Hexel& ihit = incluster[i];
      double fraction, fracTot(0.0), d2;
      for( unsigned j = 0; j < seeds.size(); ++j ) {
	fraction = 0.0;
	d2 = ( std::pow(ihit.x - centroids[j].x(),2.0) + 
	       std::pow(ihit.y - centroids[j].y(),2.0) + 
	       std::pow(ihit.z - centroids[j].z(),2.0)   )/sigma2;
	dist2[j] = d2;
	// now we set the fractions up based on hit type
	if( i == seeds[j] ) { // this cluster's seed
	  fraction = 1.0;
	} else if( isaseed[i]  ) {
	  fraction = 0.0;
	} else {	  
	  fraction = energies[j]*std::exp( -0.5*d2 );
	}
	fracTot += fraction;	
	frac[j] = fraction;
      }
      // now that we have calculated all fractions for all hits
      // assign the new fractions
      for( unsigned j = 0; j < seeds.size(); ++j ) {
	if( fracTot > minFracTot || 
	    ( i == seeds[j] && fracTot > 0.0 ) ) {
	  outclusters[i][j] = frac[j];
	} else {
	  outclusters[i][j] = 0.0;
	}	
      }
    }
    
    // save previous centroids
    prevCentroids = std::move(centroids);
    // finally update the position of the centroids from the last iteration
    centroids.resize(seeds.size());
    double diff2 = 0.0;
    for( unsigned i = 0; i < seeds.size(); ++i ) {
      centroids[i] = calculatePositionWithFraction(incluster,outclusters[i]);
      energies[i]  = calculateEnergyWithFraction(incluster,outclusters[i]);
      // calculate convergence parameters
      const double delta2 = (prevCentroids[i]-centroids[i]).perp2();
      if( delta2 > diff2 ) diff2 = delta2;
    }
    //update convergance parameter outside loop
    diff = std::sqrt(diff2);
  }
}
