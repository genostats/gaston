#ifndef _gaston_ROH_
#define _gaston_ROH_

#ifndef SHOW
#define SHOW(x) Rcpp::Rcout << #x << " = " << (x) << ", ";
#define SHOWn(x) Rcpp::Rcout << #x << " = " << (x) << std::endl;
#endif

// a class for storing the ROHS of one individual
class ROHsegments {
  public:
    // beginning and end index / positions of ROHs
    std::vector<unsigned int> begIndex;
    std::vector<unsigned int> endIndex;
    std::vector<double> begPos;
    std::vector<double> endPos;
    // nb of hets and NAs in ROHs
    std::vector<unsigned int> nbHet;

    void record(unsigned int beg, unsigned int end, double begPo, double endPo, unsigned int het) {
      begIndex.push_back(beg);
      endIndex.push_back(end);
      begPos.push_back(begPo);
      endPos.push_back(endPo);
      nbHet.push_back(het);
    }
};
 
// a class for just storing the total length of the ROH
class ROHlength {
  public:
    unsigned int nbSNPs = 0;
    unsigned int nbSegments = 0;
    double length = 0;

    void record(unsigned int beg, unsigned int end, double begPo, double endPo, unsigned int het) {
      nbSNPs += (end - beg) + 1;
      nbSegments++;
      length += (endPo - begPo);
    }
};

// ROH are defined sequences of homozygous, with a minimal number of SNPs (minNbSNPs in update)
// a minimal length (minROHLength in update)
// they can contain heterozygous SNPs, if they are are at a distance > minDistHet from the 
// extremeties of the ROH and of the others heterozygous SNPs
// NAs can be treated as homozygous or heterozygous (bool NAsAreHet)
//
// The maximal gap between SNPs is not dealt with by update(), the function going through the SNPs
// will have to call endChromosome() when such a gap is met.

template<typename T_ROH>
class ROH {
  public:
    T_ROH summary;

    bool inROH;
    unsigned int i0; // beginning index of current putative ROH
    double pos0;     // beginning position of ...
    unsigned int i1;
    double pos1;

    unsigned int lastHetIndex;  // index of last Het or NA
    double lastHetPos;          // position of last Het or NA
    unsigned int previousIndex;       // index of Homozygous SNP just before the last Het of NA
    double previousPos;         // position of Homozygous SNP just before the last Het of NA
    unsigned int het;           // nb oh heterozygous (or NAs) in current ROH

    ROH() : inROH(false), het(0) {}

    void update(unsigned int i, double pos, unsigned char gen, unsigned int minNbSNPs, double minROHLength, double minDistHet, bool NAsAreHet) {
      if(gen == 3) {
        gen = NAsAreHet?1:0;
      }
      /*
      if(debug) {
        SHOW(i);
        SHOW((int) gen);
        SHOW((int) inROH);
        SHOW(i0);
        SHOW(i1);
        // SHOW(pos0);
        // SHOW(pos1);
        SHOW(lastHetIndex);
        SHOW(lastHetPos);
        // SHOW(previousPos);
        SHOW(previousIndex);
        SHOWn(het);
      }*/
      if(gen == 1) { // heterozygote
        // pas dans ROH : rien à faire
        if(!inROH) return;
        // pas encore d'heteroz, et assez loin du bord
        // ou  bien : déjà des heteroz, et assez loin du précédent
        if( (het == 0 && pos > pos0 + minDistHet) || pos > lastHetPos + minDistHet) {
          het++;
          lastHetIndex = i;
          lastHetPos = pos;
          previousIndex = i1;
          previousPos = pos1;
          i1 = i;
          pos1 = pos;
          return;
        } 
        // Heteroz non admissible. Il faut clôturer le ROH.
        // si il n'y avait aucun heteroz, rien à faire, le ROH est juste trop court
        if(het == 0) {
          inROH = false;
          return;
        } 
        // sinon, il faut retourner à la position avant le dernier het, car il est trop près de
        // la position courante !
        i1 = previousIndex;
        pos1 = previousPos;
        het--;
        // si c'est assez long on stocke
        if( (pos1 - pos0 > minROHLength) && (i1 - i0 > minNbSNPs) ) {
          // if(debug) Rcout << "****************** push ********************\n\n";
          summary.record(i0, i1, pos0, pos1, het);
        }
        // dans tous les cas on sort du ROH courant...
        het = 0;
        inROH = false;
        return;
      } else { // homozygote
        // pas dans ROH : ouvrir ROH potentiel
        if(!inROH) {
          i0 = i1 = i;
          pos0 = pos1 = pos;
          het = 0;
          inROH = true;
          return;
        }
        // déjà dans ROH. On met à jour i1, pos1.
        i1 = i;
        pos1 = pos;
        return;
      }
    }
    // à appeler arrivé à la fin d'un chromosome ou d'un gap trop long entre les SNPs
    // la fonction clôture le ROH s'il est assez long
    void endChromosome(unsigned int minNbSNPs, double minROHLength, double minDistHet) {
      if(!inROH) return;
      // si le dernier heteroz est trop prêt du bord il faut retourner à la position précédente
      if(het > 0 && pos1 < lastHetPos + minDistHet) {
        i1 = lastHetIndex;
        pos1 = previousPos;
        het--;
      }
      // si c'est assez long on stocke
      if( (pos1 - pos0 > minROHLength) && (i1 - i0 > minNbSNPs) ) {
        // if(debug) Rcout << "****************** push ********************\n\n";
        summary.record(i0, i1, pos0, pos1, het);
      }
      // on sort du ROH !
      het = 0;
      inROH = false;
    }
};

#endif
