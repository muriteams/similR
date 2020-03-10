#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

#define a 0
#define b 1
#define c 2
#define d 3
#define row_count 4
#define col_count 5

// All statistics are stored in a map
typedef double (*funcPtr)(
    const std::vector< double > & table,
    bool normalized
);


template<typename T>
void vecprint(const std::vector< T > & x) {
  std::cout << "Here is a vector of length " << x.size() << "\n";
  for (int i = 0; i < x.size(); ++i)
    std::cout << x[i] << "\n";
  
  return;
  
}

typedef std::vector< int > vecint;

template<typename Ti, typename Tm> inline
void contingency_matrix(
    std::vector<Ti> & table,
    const Tm & M1,
    const Tm & M2,
    bool include_self,
    const vecint & exclude
  ) {
  
  // Catching error
  if ((M1.ncol() != M2.ncol()) | (M1.nrow() != M2.nrow())) 
    stop("`a` and `b` must have the same dimensions.");
  
  if (table.size() != 6)
    stop("`table` must be of size 6.");
  
  std::fill(table.begin(), table.end(), 0);
  
  int n = M1.nrow();
  int m = M1.ncol();
  
  table[row_count] = n;
  table[col_count] = m;
  
  // Set minus
  std::vector< int > rowidx(0); rowidx.reserve(n - exclude.size());
  std::vector< int > colidx(0); colidx.reserve(m - exclude.size());
  
  for (int i = 0; i < n; ++i) {
    if (std::find(exclude.begin(), exclude.end(), i) == exclude.end())
      rowidx.push_back(i);
  }
  
  for (int i = 0; i < m; ++i) {
    if (std::find(exclude.begin(), exclude.end(), i) == exclude.end())
      colidx.push_back(i);
  }

  for (auto i = rowidx.begin(); i != rowidx.end(); ++i) {
    
    for (auto j = colidx.begin(); j != colidx.end(); ++j) {
      
      if (!include_self && *i == *j)
        continue;
      
      else {
        if      ((M1(*i, *j) == M2(*i, *j)) & (M1(*i, *j) == 1)) table[a]++;
        else if ((M1(*i, *j) != M2(*i, *j)) & (M1(*i, *j) == 1)) table[b]++;
        else if ((M1(*i, *j) != M2(*i, *j)) & (M1(*i, *j) == 0)) table[c]++;
        else if ((M1(*i, *j) == M2(*i, *j)) & (M1(*i, *j) == 0)) table[d]++;
      }
      
    }
  }
}


template<typename Ti, typename Tm> inline
std::vector<Ti> contingency_matrix(
    const Tm & M1, const Tm & M2, bool include_self, const vecint & exclude
  ) {
  
  std::vector<Ti> table(6);
  contingency_matrix< Ti, Tm >(table, M1, M2, include_self, exclude);
  
  return table;
  
}

//' Contingency Table
//' @param M1,M2 Two integer matrices of the same size.
//' @param include_self Logical scalar. When `TRUE` the diagonal is
//' included in the calculation.
//' @param exclude Integer vector. List of indices to include
//' during the calculation. For example, if individual 2 needs
//' to be excluded, setting `exclude = c(2)` will include the
//' second rows and columns from calculation.
//' 
//' @export
//' 
// [[Rcpp::export(rng = false)]]
IntegerMatrix contingency_matrix(
    const IntegerMatrix & M1,
    const IntegerMatrix & M2,
    bool include_self,
    const std::vector< int > & exclude
) {
  
  if (M1.nrow() != M2.nrow())
    stop("Number of rows don't match.");
  
  if (M1.ncol() != M2.ncol())
    stop("Number of columns don't match.");
  
  std::vector< int > exclude0(exclude.size());
  for (std::vector< int >::iterator i = exclude0.begin(); i != exclude0.end(); ++i)
    (*i) -= 1;
  
  std::vector<double> table = contingency_matrix<double, IntegerMatrix>(
    M1, M2, include_self, exclude0
    );
  
  IntegerMatrix ans(2,2);
  
  ans(0, 0) = table[a];
  ans(0, 1) = table[b];
  ans(1, 0) = table[c];
  ans(1, 1) = table[d];
  
  return ans;
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Jaccard (1): `"sjaccard"` or `"jaccard"`
//' @aliases Jaccard
double sjaccard(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return table[a]/(table[a] + table[b] + table[c]);
  
}


//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Sørensen–Dice coefficient (2), Sczekanowsk (3), Nei \& Li (5): `"sdice"` or `"sczekanowsk"` or `"sneili"`
//' @aliases Sorensen-Dice
double sdice(
    const std::vector< double > & table,
    bool normalized = false
) {
  return (2.0 * table[a])/
    (2.0 * table[a] + table[b] + table[c]);
}


//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - 3w-jaccard (4): `"s3wjaccard"`
double s3wjaccard(
    const std::vector< double > & table,
    bool normalized = false
) {
  return (3.0 * table[a])/
    (3.0 * table[a] + table[b] + table[c]);
}




//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Faith (10): `"sfaith"` or `"faith"`
//' @aliases Faith
double sfaith(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return (table[a] + table[d]*0.5)/(table[a] + table[b] + table[c] + table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Gower and Legendre (11): `"sgl"` or `"gl"`
//' @aliases Gower-&-Legendre
double sgl(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return (table[a] + table[d])/(table[a] + 0.5*(table[b] + table[c]) + table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Rusell & Rao (14): `"srusrao"`
//' @aliases Rusell-&-Rao
double srusrao(
  const std::vector< double > & table,
  bool noramlized = false
) {
  return table[a]/(table[a] + table[b] + table[c] + table[d]);
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Sized Difference (24): `"dsd"` or `"sd"`
//' @aliases Sized-Difference
double dsd(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return pow(table[b] + table[c], 2.0)/
    pow(table[a] + table[b] + table[c] + table[d], 2.0);
  
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Shaped Difference (25): `"dsphd"` or `"sphd"`
//' @aliases Shape-Difference
double dsphd(
    const std::vector< double > & table,
    bool normalized = false
  ) {
  
  return (table[row_count] * (table[b] + table[c]) - pow(table[b] - table[c], 2.0))/
    pow(table[a] + table[b] + table[c] + table[d], 2.0);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Tarwid (54): `"starwid"` or `"tarwid"`.
//' @aliases Tarwid
double starwid(
    const std::vector< double > & table,
    bool normalized = false
  ) {
  
  int n = table[row_count];
  return (n*table[a] - (table[a] + table[b])*(table[a] + table[c]))/
          (n*table[a] + (table[a] + table[b])*(table[a] + table[c]));
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Pearson & Heron 1 (54): `"sph1"` or `"ph1"` or `"s14"`. This is also known as S14 in
//'    Gower and Legendre (1986).
//'    
//'    In the case of the `S14` function, following Krackhardt's 1989:
//'    
//'    \deqn{%
//'    \sqrt{\left(\frac{a}{(a + c)} - \frac{b}{(b + d)}\right)\times\left(\frac{a}{(a + b)} - \frac{c}{(c + d)}\right)}
//'    }{%
//'    S14 = [(a/(a + c) - b/(b + d))*(a/(a + b) - c/(c + d))]^(1/2)
//'    }
//'   
//'    Which is an statistic lying between 0 and 1.
//'  
//' @aliases Person-&-Heron
//' @aliases S14
double sph1(
    const std::vector< double > & table,
    bool normalized = false
  ) {

  return (table[a]*table[d] - table[b]*table[c])/
    sqrt((table[a] + table[b])*(table[a] + table[c])*(table[b] + table[d])*(table[c] + table[d]));
        
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Hamming (15): `"dhamming"` or `"hamming"`
//' @aliases Hamming
double dhamming(
    const std::vector< double > & table,
    bool normalized = false
  ) {
  
  double ans = table[b] + table[c];
  
  if (normalized) {
    return ans / (double) (table[a] + table[b] + table[c] + table[d]);
  }
    
  return ans;
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Mean Manhattan (20): `"dmh"` or `"mh"`
//'    \deqn{%
//'      D_{Mean-manhattan} = \frac{b + c}{a + b + c + d}
//'    }{%
//'     dmh = (b + c)/(a + b + c + d)
//'    }
//' @aliases Mean-Manhattan
double dmh(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return (table[b]+table[c])/
    (table[a] + table[b] + table[c] + table[d]);
    
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Dennis (44): `"sdennis"` or `"dennis"`
//' @aliases Dennis
double sdennis(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  double n = table[row_count];
  
  return (table[a]*table[d] - table[b]*table[c])/
    sqrt(n*(table[a] + table[b])*(table[a] + table[c]));
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Yuleq (61): `"syuleq"`
//' @aliases Yuleq
double syuleq(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return (table[a]*table[d] - table[b]*table[c])/(table[a]*table[d] + table[b]*table[c]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' - Yuleq (63): `"syuleqw"`
//' @aliases Yuleq
double syuleqw(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return (sqrt(table[a]*table[d]) - sqrt(table[b]*table[c]))/
    (sqrt(table[a]*table[d]) + sqrt(table[b]*table[c]));
  
}

//' @name similarity
//' @rdname similarity
//' @section Distance:
//' - Yuleq (62): `"dyuleq"`
//' @aliases Yuleq
double dyuleq(
  const std::vector< double > & table,
  bool normalized = false
) {
  
  return 2.0*table[b]*table[c]/(table[a]*table[d] + table[b]*table[c]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Michael (68):  `"smichael"` or `"michael"`
//'    
//'    \deqn{%
//'      S_{michael} = \frac{4(ad-bc)}{(a+d)^2 + (b+c)^2}
//'    }{%
//'      smichael = 4*(a*d - b*c)/[(a + d)^2 + (b + c)^2]
//'    }
//' @aliases Michael
double smichael(
  const std::vector< double > & table,
  bool normalized = false
) {
  
  return 4.0*(table[a]*table[d] - table[b]*table[c])/
    (pow(table[a] + table[d], 2.0) + pow(table[b] + table[c], 2.0));
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Dispersion (66):  `"sdisp"` or `"disp"`
//'    
//'    \deqn{%
//'      S_{Dispersion} = \frac{ad - bc}{(a + b + c + d)^2}
//'    }{%
//'      sdisp = [a*d - b*c]/(a + b + c + d)^2
//'    }
//' @aliases Dispersion
double sdisp(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return (table[a]*table[d] - table[b] * table[c])/
    pow(table[a] + table[b] + table[c] + table[d], 2.0);
  
}


//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Hamann (67):  `"shamann"` or `"hamann"`
//'    
//'    \deqn{%
//'      S_{Hamann} = \frac{(a + d) - (b + c)}{a + b + c + d}
//'    }{%
//'      shamann = [(a + d) - (b + c)](a + b + c + d)
//'    }
//' @aliases Hamann
double shamann(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  return ((table[a] + table[d]) - (table[b] - table[c]))/
    (table[a] + table[b] + table[c] + table[d]);
  
}

// Auxiliar functions for Goodman & Kruskal, and Anderberg
template<typename T>
inline double sigma(const std::vector< T > & table) {
  
  return std::max(table[a], table[b]) + std::max(table[c], table[d]) + 
    std::max(table[a], table[c]) + std::max(table[b], table[d]);
  
}

template<typename T>
inline double sigma_prime(const std::vector< T > & table) {
  
  return std::max(table[a] + table[c], table[b] + table[d]) +
    std::max(table[a] + table[b], table[c] + table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Goodman & Kruskal (69): `"sgk"` or `"gk"`
//'    
//'    \deqn{%
//'      S_{Goodman \& Kruskal} = \frac{\sigma - \sigma'}{2n - \sigma'} 
//'    }{%
//'      sgk = (s + s')/(2*n - s') 
//'    }
//'    
//'    where
//'    \eqn{\sigma = \max(a,b) + \max(c,d) + \max(a,c) + \max(b,d)}{
//'    s = max(a,b) + max(c,d) + max(a,c) + max(b,d)
//'    }, and 
//'    \eqn{\sigma' = \max(a + c, b + d) + \max(a + b, c + d)}{
//'    s' = max(a + c, b + d) + max(a + b, c + d)
//'    }
//'   
//' @aliases Goodman-&-Kruskal
double sgk(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  double s = sigma(table);
  double s_prime = sigma_prime(table);
  
  return (s - s_prime)/(2.0 * ((double) table[row_count]) - s_prime);
  
}

//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Anderberg (70): `"sanderberg"` or `"anderberg"`
//'    
//'    \deqn{%
//'      S_{Anderberg} = \frac{\sigma - \sigma'}{2n} 
//'    }{%
//'      sanderberg = (s + s')/(2*n) 
//'    }
//'    
//'    where \eqn{\sigma}{s} and \eqn{\sigma}{s'} are defined as in (69).
//'   
//' @aliases Anderberg
double sanderberg(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  double s = sigma(table);
  double s_prime = sigma_prime(table);
  
  return (s - s_prime)/(2.0 * ((double) table[row_count]));
  
}





//' @name similarity
//' @rdname similarity
//' @section Similarity:
//' -  Peirce (73): `"speirce"` or `"peirce"`
//'    
//'    \deqn{%
//'      S_{Peirce} = \frac{ab + bc}{ab + 2bc + cd}
//'    }{%
//'      speirce = (a*b + b*c)/(a*b + 2*b*c + c*d)
//'    }
//'   
//' @aliases Peirce
double speirce(
  const std::vector< double > & table,
  bool normalized = false
) {
  
  return (table[a]*table[b] + table[b]*table[c])/
    (table[a]*table[b] + 2*table[b]*table[c] + table[c]*table[d]);
  
}

//' @name similarity
//' @rdname similarity
//' @aliases Fscore
//' @section Similarity: 
//' In the case of `fscore`, ask Kyosuke Tanaka.
double fscore(
    const std::vector< double > & table,
    bool normalized = false
) {
  
  double precision = table[a]/(table[a] + table[c]);
  double recall    = table[a]/(table[b] + table[a]);
  
  return 2.0 * (precision * recall / (precision + recall));
  
}

// -----------------------------------------------------------------------------

typedef double (*funcPtr)(
    const std::vector< double > & table,
    bool normalized
  );

// [[Rcpp::export(rng = false)]]
IntegerMatrix reduce_dim(IntegerMatrix & x, int k) {
  
  // No reduction, return the same
  if (k == -1)
    return x;
  
  // Preparing reduced version of the matrix
  IntegerMatrix ans(x.nrow() - 1, x.ncol() - 1);
  ans.fill(0);
  
  int i0=0, j0;
  for (int i = 0; i < x.nrow(); ++i) {
    
    if (i == k)
      continue;

    j0 = 0;    
    for (int j = 0; j < x.ncol(); ++j) {
      
      // Self
      if (j == k) 
        continue;
      
      if (x(i,j) == 1)
        ans(i0, j0) = 1;

      // Out of range (can't do more)      
      if (++j0 == ans.ncol())
        break;
     
    }
    
    // Out of range (can't do more)
    if (++i0 == ans.nrow())
      break;
  }
  
  return ans;
  
}

void getmetric(std::string s, funcPtr & fun) {
  
  if      ((s == "sjaccard") | (s == "jaccard"))      fun = &sjaccard;
  else if ((s == "sdice") | (s == "sczekanowsk") | (s == "sneili")) fun = &sdice;
  else if ((s == "s3wjaccard") | (s == "3wjaccard"))  fun = &s3wjaccard;
  else if ((s == "sfaith") | (s == "faith"))          fun = &sfaith;
  else if ((s == "sgl") | (s == "gl"))                fun = &sgl;
  else if ((s == "srusrao") | (s == "rusrao"))        fun = &srusrao;
  else if ((s == "sph1") | (s == "ph1") | (s == "s14")) fun = &sph1;
  else if ((s == "dhamming") | (s == "hamming"))      fun = &dhamming;
  else if ((s == "dennis") | (s == "sdennis"))        fun = &sdennis;
  else if ((s == "starwid") | (s == "tarwid"))        fun = &starwid;
  else if (s == "syuleq")                         fun = &syuleq;
  else if (s == "syuleqw")                        fun = &syuleqw;
  else if (s == "dyuleq")                         fun = &dyuleq;
  else if ((s == "smichael") | (s == "michael"))      fun = &smichael;
  else if ((s == "speirce") | (s == "peirce"))        fun = &speirce;
  else if ((s == "sgk") | (s == "gyk"))               fun = &sgk;
  else if ((s == "sanderberg") | (s == "anderberg"))  fun = &sanderberg;
  else if ((s == "shamann") | (s == "hamann"))        fun = &shamann;
  else if ((s == "dmh") | (s == "mh"))                fun = &dmh;
  else if ((s == "dsd") | (s == "sd"))                fun = &dsd;
  else if ((s == "dsphd") | (s == "sphd"))            fun = &dsphd;
  else if ((s == "sdisp") | (s == "disp"))            fun = &sdisp;
  else if ((s == "fscore") | (s == "Fscore"))         fun = &fscore;
  else Rcpp::stop("The statistic '%s' is not defined.", s);
  
  return ;
}

// Applies whatever similarity/distance metric should be applying to all the
// requested combinations
// [[Rcpp::export(name=".similarity", rng = false)]]
NumericMatrix similarity(
    const ListOf<IntegerMatrix> & M,
    const std::vector< std::string > & statistic,
    bool normalized   = false,
    bool firstonly    = false,
    bool include_self = false,
    bool exclude_j    = false
) {
  
  int N  = M.size();
  int NN = firstonly? N - 1: N * (N - 1) / 2;
  int nfuns = statistic.size();
  NumericMatrix ans(NN, nfuns);
  NumericVector I(NN), J(NN);
  
  int pos = 0;
  
  std::vector< int > exclude(0u);
  
  int firstloop = firstonly ? 1 : N;  
  
  if (exclude_j)
    exclude.push_back(0u);
  
  std::vector< double > ctable(6u);
  
  for (int i = 0; i < firstloop; ++i) {
    for (int j = i; j < N; ++j) {
      
      if (exclude_j)
        exclude[0] = j;
      
      if (i == j)
        continue;
      int s = 0;
      I[pos] = i + 1;
      J[pos] = j + 1;
      
      // Computing contingency table, this will be used for all measurements.
      contingency_matrix(ctable, M[i], M[j], include_self, exclude);
      
      for (auto fname = statistic.begin(); fname != statistic.end(); ++fname) {
        
        // Getting the function
        funcPtr fun;
        getmetric(*fname, fun);
      
        ans(pos, s++) = fun(ctable, normalized);
      }
      pos++;
    }
  }
  
  
  return cbind(I,J,ans);
  
}

// No longer needed
#undef a
#undef b
#undef c
#undef d
#undef row_count
#undef col_count

