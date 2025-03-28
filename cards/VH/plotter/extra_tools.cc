//////////////////////////////////////////////////////////////////////
//	Author: Mohammad Abrar Wadud, Univeristy of Minnesota			//
//	Date: Jan/03/2020												//
//////////////////////////////////////////////////////////////////////

#ifndef EXTRATOOLS_H
#define EXTRATOOLS_H

// #include "KDEProducer1D.h"
#include "math.h"
#include "dirent.h"
#include "errno.h"
#include "dirent.h"
#include "stdio.h"
#include "cassert"
#include "map"
#include "cmath"
#include "sys/stat.h"
#include "Math/DistFunc.h"
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/replace.hpp"
#include "boost/algorithm/string_regex.hpp"
#include "boost/type_index.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/lexical_cast/bad_lexical_cast.hpp"
#include "iostream"
#include "sstream"
#include "fstream"
#include "vector"
#include "algorithm"
#include "cctype"
#include "locale"
#include "numeric"
#include "regex"
#include "iterator"
#include "string"
#include "chrono"
#include "ctime"
#include "limits"
#include "cstdio"
#include "memory"
#include "stdexcept"
#include "array"
#include "iomanip"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include "TKey.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TAxis.h"
#include "TObjArray.h"
#include "TColor.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TCollection.h"
#include "TKey.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMarker.h"
#include "TColor.h"
#include "TLatex.h"
#include "RtypesCore.h"
#include "TEventList.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"

#ifndef sanityCheck
#define sanityCheck(BoolCondition) do{ if(!(BoolCondition)){ std::cout<<"\nError! Bug in file "<<__FILE__<<"\tline "<<__LINE__<<"\nAborting..."<<std::endl; abort();}} while (0)
#endif

#ifndef printLine
#define printLine(){std::cout<<__FILE__<<"\tline "<<__LINE__<<std::endl;}
#endif

/*************************************************************Declarations*************************************************************/
Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t deltaPhi(Double_t phi1, Double_t phi2);
template <typename anytype1, typename anytype2>
void setBit(anytype1 & _container, anytype2 _pos, Bool_t _bitVal);
template <typename anytype1, typename anytype2>
Bool_t getBit(anytype1 _container, anytype2 _pos);
std::string removeNonAlpha(std::string word); /// Removes non-alphanumeric characters from a string
template <class any_number>
std::string removeTrailingZeros(any_number number); //// Converts 1.12000000 to a string "1.12"
Bool_t file_exists(std::string fileName);
Int_t mkdir(std::string dir_path);
std::vector<std::string> getObjectList(std::string filepath, std::string objtype, std::vector<std::string> exclusion_list= {});
Bool_t match(std::string _needle, std::string _haystack); //// Checks if a wildcard (e.g. "*a*b*") exists in "abrar")
Bool_t matchManyNeedles(std::vector<std::string> _needles, std::string _haystack);
std::string ReadNthLine(std::string filename, int N);
UInt_t countLines(std::string filename);
std::vector<std::string> split_string(std::string _string, std::string _delimiter=",", Bool_t _trim=1);
std::string get_cell(std::string filename, UInt_t row, UInt_t column, std::string _delimiter=",");
std::string getFileName(std::string _filepath);
std::string getDirPath(std::string _somePath);
std::vector<std::string> getNonemptyLines(std::string filepath);
std::vector<std::string> getNonemptyLinesWithFilenameKeyword(std::string filepath, std::string keyword, std::string exclude="");
std::vector<std::string> getLinesRegex(std::string _filepath, std::string _regexStr);
TH1* getHistFromFile(std::string _histName, std::string _filename, Bool_t _verbose=1, TDirectory *_dir = nullptr);
TObject *getObjectFromFile(std::string _objectName, std::string _filename);
void removeBinErrors(TH1* _hist);
TH1* rebinNHist(TH1* _hist, Int_t _N);
TH1* rebinHist(TH1* _hist, std::vector<Double_t> _newBins);
TH1* rebinHist(TH1* _hist, TH1* _templateHist);
std::vector<Double_t> getGoodBins(TH1* _hist, Double_t _statUnc, Double_t _reScale = -999., Int_t _nbinPar=10);
Double_t sumNextNbins(TH1* _hist, Int_t _n, Int_t _curr);
void copyHistAtts(TH1* _source, TH1* _mock);
Double_t getSumW(std::string _cutflowfile);
void ltrim(std::string &s);
void rtrim(std::string &s);
void trim(std::string &s);
std::string ltrim_copy(std::string s);
std::string rtrim_copy(std::string s);
std::string trim_copy(std::string s);
Double_t ams(Double_t _s, Double_t _b);
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 0);
std::string getUnit(TH1* _hist);
std::string getUnit(std::string ystring);
std::string eraseUnit(std::string ystring);
std::string first_numberstring(std::string const & str);
std::vector<Double_t> getXbins(TH1 *hist);
void closeTChain(TChain * & _chain);
void setFrameColor(TAxis* _axis, std::string _color);
void setFrameColor(TH1* _hist, std::string _color);
void setFrameColor(THStack* _stack, std::string _color);
TTree *loadTreeFromFile(std::string _treeName, std::string _filename);
Char_t isDirectory(std::string filePath, Bool_t _verbose=0);
TGraph* removeErrors(TGraphErrors* graphToConvert);
void addPointToGraph(TGraphAsymmErrors * _graph, Double_t _x, Double_t _xErrL, Double_t _xErrH, Double_t _y, Double_t _yErrL, Double_t _yErrH);
void addPointToGraph(TGraph * _graph, Double_t _x, Double_t _y);
void addPointToGraph(TGraphErrors * _graph, Double_t _x, Double_t _ex, Double_t _y, Double_t _ey);
void graphStats(TGraphAsymmErrors* graph, Double_t &mean, Double_t &stdev);
TH1D* graph2hist(TGraphAsymmErrors* graph, UInt_t ndivs, Double_t ylow, Double_t yhigh);
Int_t writeToFile(TObject *_object, std::string _filePath, std::string _mode = "RECREATE", Bool_t _verbose=1, std::string _writeName = "", std::string _tdir = "");
Int_t writeHistToFile(TH1 *_object, std::string _filePath, std::string _mode = "RECREATE", Bool_t _verbose=1, std::string _writeName = "");
std::string findAndReplaceAll(std::string data, std::string toSearch, std::string replaceStr);
TChain *openTree(std::string _inFile, std::string _treeName="", Bool_t _verbose = 1);
TChain *openTChain(std::string _chainListFile, std::string _treeName="", Bool_t _verbose = 1);
TChain *openTChain(std::vector<std::string> _chainList, std::string _treeName="", Bool_t _verbose = 1);
// TChain *openTChainWithFilesInDir(std::string _dirPath, std::string _treeName="");
std::vector<std::pair<std::string, std::string>> getBranchList(std::string _treePath, std::string _treeName="", Bool_t _verbose=1);
std::vector<std::string> listFilesInDir(std::string _dirPath, std::string _regexStr="", Bool_t _verb=0);
Bool_t matchRegex(std::string _str,std::string _regexStr);
std::string getTreeNameInFile(std::string _filePath);
TH1F *mergeBins(std::string _fileList, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol=0, Int_t _xSecCol=2, std::string _path="");
TH1F *mergeBins(std::vector<std::string> _fileList, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol=0, Int_t _xSecCol=2, std::string _path="");
std::string vLookup(std::string _lookupKey, std::string _inFile, Int_t _lookupCol, Int_t _valCol, Bool_t _regex=0, std::string _delimiter=",", Bool_t _silent = 1);
std::string hvLookup(std::string _rowLookupKey, std::string _colLookupKey, std::string _inFile, Int_t _lookupRow=0, Int_t _lookupCol=0, Bool_t _regex=0, std::string _delimiter = ",", Bool_t _silent = 1);
Bool_t isROOTfile(std::string _filepath);
std::vector<Float_t> getXlimits(std::vector<TH1*> _hists, Float_t _binThreshold=0.);
void clearHeap();
Double_t weightedYmean(TH1 *_hist);
Bool_t branchExists(std::string _branchName, TTree *_tree);
Float_t getMean(std::vector<Float_t> _set);
template <class ObjType>
ObjType copyObjectDeleteSrc(ObjType *_original);
template<typename T1, typename T2>
Int_t findSecondaryIndex(T1 searchIndex, std::vector<T2> container);
Short_t findSecondaryIndex(Short_t searchIndex, std::vector<Short_t> container);
std::string getCurrentTime();
template <typename anytype>
void eraseElement(std::vector<anytype> & _vec, UInt_t _Index2Erase);
Double_t getCategoryBoundary(TH1*_signal, TH1*_background);
std::vector<Double_t> vecString2vecDouble(std::vector<std::string> _numStrings);
Bool_t stringIsInteger(std::string _isThisAnInt);
Bool_t stringIsNumber(std::string _isThisANumber);
TDirectory *mkTFileDir(TFile *_file, std::string _dir);
Int_t caselessSubstringPos( std::string str1, std::string str2);
Bool_t findStringIC(const std::string & strHaystack, const std::string & strNeedle);
Bool_t stringIsEmpty(std::string str);
template<typename T>
UChar_t setLineAtts(T* graph, std::string _atts, std::string _delimiter=",");
template<typename T>
UChar_t setMarkerAtts(T* graph, std::string _atts, std::string _delimiter=",");
template<typename T>
UChar_t setFillAtts(T* graph, std::string _atts, std::string _delimiter=",");
UChar_t setAxisAtts(TAxis* _axis, std::string _atts, std::string _delimiter=",");
template<typename T>
UChar_t setPadAtts(T* _pad, std::string _atts, std::string _delimiter=",");
template<typename T>
void setFrameAtts(T* _hist, std::string _style, std::string _del1=";", std::string _del2=",");
template<typename T>
void setHistStyle(T* _hist, std::string _style, std::string _del1=";", std::string _del2=",");
std::string stringToUpper(std::string _str);
std::string stringToLower(std::string _str);
template<typename T1, typename T2>
Int_t findIndex(const std::vector<T1> & _haystack, T2 _needle);
Bool_t objectExists(TFile *_File, std::string _objectName);
Bool_t objectExists(std::string _FilePath, std::string _objectName);
std::vector<Float_t> strToFloatList(std::string _listString, std::string _delimiter=",");
std::vector<Double_t> strToDoubleList(std::string _listString, std::string _delimiter=",");
template<typename T>
void setAxes(T* _graph, std::string _axesOptions) {};
Int_t hex2rootColor(std::string _hexCode);
std::string sysexec(std::string cmd);
std::vector<Double_t> getNpercentMinInterval(TH1 *_hist, Float_t _N);
std::vector<Double_t> getNpercentLinearInterval(TH1 *_hist, Float_t _N);
std::vector<std::string> prefixVecString(std::vector<std::string> _vecStr, std::string _prefix, Int_t _insPos = 0);
void normalizeHist(TH1* _hist, Double_t _norm=100., Bool_t _full = 0);
void normalizeHist(TH1& _hist, Double_t _norm=100., Bool_t _full = 0);
std::vector<Int_t> strToIntList(std::string _listString, std::string _delimiter=",");
Double_t spearmanR(TH2 *_hist, Int_t _Nbins = 1000);
Float_t relInvMass(Float_t _pT1, Float_t _eta1, Float_t _phi1, Float_t _pT2, Float_t _eta2, Float_t _phi2);
Float_t relTransMass(Float_t _pT1, Float_t _phi1, Float_t _pT2, Float_t _phi2);
template<typename T>
void remove_intersection(std::vector<T>& a, std::vector<T>& b);
TH1* getSqrtHist(TH1* _hist);
Bool_t hasNegativeBins(TH1* _hist, Float_t _thres = 0.);
void removeNegativeBins(TH1* _hist, Float_t _thres = 0.);
TH1* copyHistSubRange(TH1* _hist, Float_t _xMin, Float_t _xMax, std::string newName = "");
Float_t foldedPhi(Float_t _phi);
bool isInteger(const std::string & s);
Char_t RGB2ColPlt(std::string _CSVFile, Int_t _firstCol);
TH1* getAbsHist(TH1* _hist);
template <typename T>
std::string ToSciString(T value, UShort_t nDec = 2, std::string powFormat = "e");
template <typename T>
std::string ToSciString(T value, T valErr, UShort_t nDec, std::string pmSign = "#pm", std::string powFormat="e");
void cmsLegend(TPad* _pad, std::string _runInfo, std::string _2ndary="", Float_t _scale=1.);
void TGraphErrorsToCSV(TGraphErrors *_graph, std::string _CSVoutfile);
void clearTHStack(THStack *hStack);
Float_t getPunziSmin(Float_t _B, Float_t _a, Float_t _b);
TH1* getPunziSminHist(TH1* _Bhist, Float_t _a, Float_t _b);
TGraphErrors* makeROC(TH1* _sigEffHist, TH1* _bgEffHist);
TH1* getMaxHist(std::vector<TH1*> hists, Bool_t onlyAbs = 1);
// TH1* getInverseHist(TH1* _hist);
void	drawMaxContour(TH2* _hist, Double_t _contour, std::string _contCol="FF0000", Float_t _textSize=0.027, Float_t _contLineWidth=4, Float_t _markerSize = 2.);
void	drawMinContour(TH2* _hist, Double_t _contour, std::string _contCol="FF0000", Float_t _textSize=0.027, Float_t _contLineWidth=4, Float_t _markerSize = 2.);
Double_t separationTMVA(const TH1* _histSig, const TH1* _histBG, Bool_t _normalize = 1, Bool_t _preScaledByWidth = 0);
Float_t getMaxInRange(TH1* _hist, Float_t _xMin, Float_t _xMax);
bool TestFitSuccess(bool verbose = false);
std::string randStr(UInt_t _len=6);
std::string randHexStr(UInt_t _len=6);
TGraph* subtractTGraphs(TGraph* g1, TGraph* g2);
TGraphErrors* subtractTGraphErrors(TGraphErrors* g1, TGraphErrors* g2, Bool_t g1err, Bool_t g2err);
TGraph* divideTGraphs(TGraph* g1, TGraph* g2);
TGraphErrors* divideTGraphErrors(TGraphErrors* g1, TGraphErrors* g2, Bool_t g1err, Bool_t g2err);
void divideWithoutError(TH1* _mainHist, TH1* _divisor);
void divideSelf(TH1* _mainHist);
std::pair<Double_t, Double_t>  				getHistQuantile(TH1* _hist, Double_t _quantile);
std::pair<Double_t, Double_t>  				getHistQuantile2D(TH2* _hist, Double_t _quantile, Int_t _bin, Int_t _axis = 0);
void 	addContent(TH1D* _mainHist, TH1D* _addErrorsHist);
TH1D* 	getEnvelope(std::vector<TH1D*> _systHists, TH1D* _nominalHist, std::string _newName);
TH1D* 	getSymmetricEnvelope(std::vector<TH1D*> _systHists, std::string _newName);
TH1D* 	addError(TH1D* _mainHist, TH1D* _addErrorsHist, std::string _newName);
void 	addStatErrorToContent(TH1D* _variationHist, TH1D* _statErrorsHist, Float_t _addSign = 1., Bool_t _truncateNegative = 1) ;
std::vector<TH1D*>	getPDFHessianSystSum(std::vector<TH1D*> _systHists, TH1D* _nominalHist, std::string _newName, Bool_t _addAlphaS, Float_t _rescaleAlphaSunc);
TH1D* 	getPDFHessianSyst(std::vector<TH1D*> _systHists, TH1D* _nominalHist, std::string _newName);
std::vector<TH1D*>												getUpDnEnvelope(std::vector<TH1D*> _systHists, std::string _newName);
TH1D* 	addEnvelope(std::vector<TH1D*> _systHists, TH1D* _nominalHist, TH1D* _prevEnvelop, std::string _newName, Bool_t _addQuad = 0);
template <class HISTCLASS>
std::vector<Double_t> getHistMinMax(std::vector<HISTCLASS*>  _hist, Bool_t _ignoreZero = 1, Bool_t _ignoreError = 0);
void scaleHistErrBandByTanh(TH1* _hist, Double_t _yCentre);
TGaxis* makeTanhYaxis(TPad* _pad, Double_t _yMin, Double_t _yMax, Double_t _yCentre);
void moveOverflowToLastBin(TH1* _hist);
void moveUnderflowToFirstBin(TH1* _hist);
Double_t getWeightedBinVariation(std::vector<TH1*> _variations, TH1* _nominalHist);
Double_t getHistSupremum(TH1* _hist, const Double_t _upperBound = std::numeric_limits<Double_t>::max()/1000000.);
Double_t getHistInfimum(TH1* _hist, const Double_t _lowerBound = 1000000.*std::numeric_limits<Double_t>::min());
void deleteObjectFromTFile(std::string _filePath, std::string _objName);
void copyFile(std::string _sourcePath, std::string _destPath);
void writeStringToFile(std::string _stringToWrite, std::string _filePath, Bool_t _append = 0, Bool_t _verbose = 0);
TGraphAsymmErrors* getPoissonErrorRatioGraph(TH1* _dataHist, TH1* _predHist, Float_t _confidence = 0.682689492137086);
TGraphAsymmErrors* getPoissonErrorGraph(TH1* _hist, Bool_t _divideByWidth=0, Float_t _confidence = 0.682689492137086);
template <typename T>
Bool_t areFloatsEqual(T val1, T val2, T tmpEpsilon);
template <typename T>
Bool_t areVectorsEqual(std::vector<T> vec1, std::vector<T> vec2, T tmpEpsilon);
TGraphAsymmErrors* systematicBandFromHists(TH1* nomHist, TH1* upHist, TH1* dnHist, Bool_t skipZero = 1);

TGraphAsymmErrors* systematicBandFromHists(TH1* nomHist, TH1* upHist, TH1* dnHist, Bool_t skipZero){

	TGraphAsymmErrors* systBandGraph = new TGraphAsymmErrors();

	for(Int_t i = 1; i <= nomHist->GetNbinsX(); i++){

		Double_t yNominal = nomHist->GetBinContent(i);

		if(skipZero && (yNominal == 0.)) continue;
		
		Double_t yMax = std::max(upHist->GetBinContent(i), dnHist->GetBinContent(i));
		Double_t yMin = std::min(upHist->GetBinContent(i), dnHist->GetBinContent(i));

		addPointToGraph(systBandGraph, nomHist->GetXaxis()->GetBinCenter(i), nomHist->GetXaxis()->GetBinWidth(i)/2., nomHist->GetXaxis()->GetBinWidth(i)/2.,
			yNominal, yMax - yNominal, yNominal - yMin);
	}

	return systBandGraph;
};
// template<typename T1, typename T2>
// std::vector<T2> getMapValueVec(const std::map<T1, T2> & theMap);


// void setTFileDir(std::vector<TH1*> _hists);

// void setTFileDir(std::vector<TH1&> _hists, TFile & _file, std::string _dir = ""){
// 	_file.mkdir(_dir.c_str());
// 	for(TH1 & iHist : _hists){
// 		iHist.SetDirectory(_file.GetDirectory(_dir.c_str()));
// 	}
// };

// std::string getTFormulaParam(TF1* _func){

// };


template <typename anytype>
struct TTreeReaderAnyValue;
template <typename anytype>
struct TTreeReaderVectorValue;
template <typename anytype>
struct TTreeReaderArrayValue;
struct TTreePathok;
struct plot_variable;
struct histogram_template;
struct twoDhistogram_template;
struct bit_histogram_template;
struct bit_twoDhistogram_template;
template <typename anytype>
struct vector_association;
struct BinCollection;
struct signal_atts;
struct sample;
struct Profile2D;
struct parseOptions;
struct bitBar2D;
struct kFactorMap;
struct effectiveAreaMap;
struct isoPtScalingMap;
struct coronaCorrections;
struct BGset;
template <typename anytype>
struct smarter_ptr;
template <typename anytype>
struct indexer;
struct DinkyHister;
struct logStream;
struct correlationMatix;
struct scaleFactor1D;
struct pTEtaScaleFactor2D;
struct etaPtScaleFactor2D;
struct binCenterTracker;
struct QuantileCorrector;




class CSVReader;
class PileupReWeighting;


const std::map<std::string, Double_t> xsec_unit_map = {
	{"fb", 1.0e-3},
	{"pb", 1.0},
	{"nb", 1.0e3}
};


namespace BayesianBlocks {
////// Adapted from https://github.com/gipert/bayesian-blocks/tree/master/cpp

// handy aliases
	namespace bb {

// data containers
		using array         = std::vector<double>;
		using data_array    = std::vector<double>;
		using weights_array = std::vector<double>;
		using pair          = std::pair<double, double>;

// time
		using clock = std::chrono::high_resolution_clock;
		using std::chrono::duration_cast;
		using us = std::chrono::microseconds;
	}

// core utility
	bb::array blocks(bb::data_array data, bb::weights_array weights, const double p = 0.01,
		bool counter = false, bool benchmark = false);

	bb::array blocks(bb::data_array data, const double p = 0.01,
		bool counter = false, bool benchmark = false);

// rebin a ROOT histogram
	TH1* rebin(TH1* h_in, const double p = 0.01,
		bool counter = false, bool benchmark = false);

// core utility
	bb::array blocks(bb::data_array data, bb::weights_array weights, const double p, bool counter, bool benchmark) {

		auto start = bb::clock::now();

	// sanity checks
		if (data.size() != weights.size()) {
			throw std::domain_error("ERROR: data and weights vectors are of different sizes");
		}

		if (data.size() == 0) {
			throw std::invalid_argument("ERROR: empty arrays provided as input");
		}


		if (std::unique(data.begin(), data.end()) != data.end()) {
			throw std::invalid_argument("ERROR: duplicated values found in input");
		}

		const auto N = data.size();

	// sort and copy data
		std::vector<bb::pair> hist; hist.reserve(N);
		for (std::size_t i = 0; i < N; ++i) hist.emplace_back(data[i], weights[i]);
			std::sort(hist.begin(), hist.end(), [](bb::pair a, bb::pair b) { return a.first < b.first; });
		for (size_t i = 0; i < N; ++i) {
			data[i]    = hist[i].first;
			weights[i] = hist[i].second;
		}

	// build up array with all possible bin edges
		bb::array edges(N+1);
		edges[0]   = data[0];
		for (std::size_t i = 0; i < N-1; ++i) edges[i+1] = (data[i]+data[i+1])/2.;
			edges[N]   = data[N-1];

		assert(std::unique(edges.begin(), edges.end()) == edges.end());

	// let's use here Cash statistics and calibrated prior on number of change points
		auto cash = [](double N_k, double T_k) { return N_k * std::log(N_k/T_k); };
		auto ncp_prior = std::log(73.53 * p * std::pow(N, -0.478)) - 4;

	// arrays to store results
		bb::array last(N);
		bb::array best(N);

		auto init_time = bb::duration_cast<bb::us>(bb::clock::now() - start).count();
		start = bb::clock::now();

	// do the actual recursive computation
		for (std::size_t k = 0; k < N; ++k ) {
			bb::array A(k+1);
			for (std::size_t r = 0; r == 0 or r <= k; ++r) {
				A[r] = cash(
					std::accumulate(weights.begin()+r, weights.begin()+k+1, 0),
					edges[k+1] - edges[r]
					) + ncp_prior + (r == 0 ? 0 : best[r-1]);
			}
			last[k] = std::distance(A.begin(), std::max_element(A.begin(), A.end()));
			best[k] = *(std::max_element(A.begin(), A.end()));

			if (counter) std::cout << '\r' << k << '/' << N << std::flush;
		}
		if (counter) std::cout << std::endl;

		auto loop_time = bb::duration_cast<bb::us>(bb::clock::now() - start).count();
		start = bb::clock::now();

		std::vector<int> cp;
		for (auto i = N; i != 0; i = last[i-1]) cp.push_back(i); cp.push_back(0);

			std::reverse(cp.begin(), cp.end());
		bb::array result(cp.size(), 0);
		std::transform(cp.begin(), cp.end(), result.begin(),
			[edges](size_t pos) { return edges[pos]; });

		auto end_time = bb::duration_cast<bb::us>(bb::clock::now() - start).count();
		if (benchmark) {
			std::cout << "init: ";
			init_time > 1000 ?
			std::cout << init_time/1.E3 << " s" :
			std::cout << init_time      << " us";
			std::cout << std::endl;
			std::cout << "loop: ";
			loop_time > 1000 ?
			std::cout << loop_time/1.E3 << " s" :
			std::cout << loop_time      << " us";
			std::cout << std::endl;
			std::cout << "end: ";
			end_time > 1000 ?
			std::cout << end_time/1.E3 << " s" :
			std::cout << end_time      << " us";
			std::cout << std::endl;
		}
		return result;
	}

	bb::array blocks(bb::data_array data, const double p,
		bool counter, bool benchmark) {
		std::map<double, int> hist;
		for (auto& i : data) {
			if (hist.find(i) == hist.end()) hist.emplace(i, 1);
			else hist[i]++;
		}
		bb::data_array x;
		bb::weights_array weights;
		for (auto& i : hist) {
			x      .push_back(i.first);
			weights.push_back(i.second);
		}
		return BayesianBlocks::blocks(x, weights, p, counter, benchmark);
	}

	TH1* rebin(TH1* h_in, const double p, bool counter, bool benchmark) {
		const auto Nb = h_in->GetNbinsX();
		bb::data_array x;
		bb::weights_array weights;
		bb::array edges;
		int i_first = 1;
		for (int i = 1; i < Nb; ++i ) {
			if (h_in->GetBinContent(i) >0) {
				edges  .push_back(h_in->GetBinLowEdge(i));
				x      .push_back(h_in->GetBinCenter(i));
				weights.push_back(h_in->GetBinContent(i));
				i_first = i;
				break;
			}
		}
		for (int i = i_first+1; i < Nb; ++i ) {
			auto c = h_in->GetBinContent(i);
			if (c == 0) continue;
			edges  .push_back(h_in->GetBinLowEdge(i));
			x      .push_back(h_in->GetBinCenter(i));
			weights.push_back(c);
		}
		edges  .push_back(h_in->GetBinLowEdge(Nb));
		x      .push_back(h_in->GetBinCenter(Nb));
		weights.push_back(h_in->GetBinContent(Nb));
		edges  .push_back(h_in->GetXaxis()->GetBinUpEdge(Nb));

		auto result = BayesianBlocks::blocks(x, weights, p, counter, benchmark);
		bb::array reAlignedBins;

		std::sort(edges.begin(), edges.end());
		for (UInt_t i = 0 ; i < result.size(); i++) {
			std::vector<double>::iterator iNearestBinUpEdge = std::upper_bound(edges.begin(), edges.end(), result[i]);
			if (i == 0) {
				std::vector<double>::iterator firstBinLowEdge = iNearestBinUpEdge;
				firstBinLowEdge--;
				reAlignedBins.push_back(*firstBinLowEdge);
			}
			reAlignedBins.push_back(*iNearestBinUpEdge);
		}

		auto h_out = new TH1D((std::string(h_in->GetName()) + "_b").c_str(), h_in->GetTitle(), reAlignedBins.size()-1, reAlignedBins.data());
		for (int b = 1; b < h_in->GetNbinsX(); ++b) {
			auto c = h_in->GetBinContent(b);
			auto bin = h_out->FindBin(h_in->GetBinCenter(b));
			h_out->SetBinContent(bin, h_out->GetBinContent(bin) + c);
		}
	// h_out->SetBinContent(0, h_in->GetBinContent(0));
	// h_out->SetBinContent(result.size()-1, h_in->GetBinContent(h_in->GetNbinsX()+1));
	// h_out->Scale(1, "width");
		return h_out;
	}
}


/*************************************************************Definitions*************************************************************/

template <typename anytype1, typename anytype2>
void setBit(anytype1 & _container, anytype2 _pos, Bool_t _bitVal) {
	_container ^= (-_bitVal ^ _container) & (1UL << _pos);
};


template <typename anytype1, typename anytype2>
Bool_t getBit(anytype1 _container, anytype2 _pos) {
	return (_container>>_pos) & 1;
};


template <typename anytype>
struct TTreeReaderAnyValue {
	TTreeReaderValue<anytype> *val = nullptr;
	TTreeReaderAnyValue(TTreeReader & ttreereader, std::string branchname, Bool_t _verbose = 0) {
		set(ttreereader, branchname, _verbose);
	};
	TTreeReaderAnyValue() {};
	~TTreeReaderAnyValue() {
		delete val;
		val = nullptr;
	};
	void set(TTreeReader & ttreereader, std::string branchname, Bool_t _verbose = 0) {

		// if(!branchExists(branchname, ttreereader.GetTree())){
		// 	std::cout<<"\t\tError! Branch "<<branchname<<" does not exist in TTree "<<ttreereader.GetTree()->GetName()<<std::endl;
		// 	exit(EXIT_FAILURE);
		// }
		
		delete val;

		// if (!branchExists(branchname, ttreereader.GetTree())) {
		// 	std::cerr << "Error! Branch '" << branchname << "' not found!" << '\n';
		// 	exit(EXIT_FAILURE);
		// }

		val = new TTreeReaderValue<anytype>(ttreereader,branchname.c_str());

		if (_verbose) std::cout<<"\t\tInitialized branch <"<< boost::typeindex::type_id<anytype>().pretty_name() <<"> "<<branchname<<std::endl;
	};
	operator anytype () const {

		// if (val->GetSetupStatus() < 0) {
		// 	std::cerr << "Error " << val->GetSetupStatus() << " setting up TTreeReader branch for " << val->GetBranchName() << '\n';
		// 	exit(EXIT_FAILURE);
		// }

		// anytype iVal = **val;

		return **val;
	}

	anytype get() {

		return **val;
	};
};


template <typename anytype>
struct TTreeReaderVectorValue: TTreeReaderAnyValue<std::vector<anytype>> {
	TTreeReaderVectorValue(TTreeReader & ttreereader, std::string branchname, Bool_t _verbose = 1):TTreeReaderAnyValue<std::vector<anytype>>(ttreereader, branchname, _verbose) {
	};
	TTreeReaderVectorValue() {};
	~TTreeReaderVectorValue() {
	};
	anytype operator[](UInt_t index) {
		anytype iVal = (*(this->val))->at(index);
		return iVal;
	};

	anytype at(UInt_t index) {
		anytype iVal = (*(this->val))->at(index);
		return iVal;
	};
	UInt_t size() {
		return (*(this->val))->size();
	};

	//define begin & end function for iterator
};


template <typename anytype>
struct TTreeReaderArrayValue {
	TTreeReaderArray<anytype> *val = nullptr;
	TTreeReaderArrayValue(TTreeReader & ttreereader, std::string branchname, Bool_t _verbose = 1) {
		set(ttreereader, branchname, _verbose);
	};
	TTreeReaderArrayValue() {};
	~TTreeReaderArrayValue() {
		delete val;
		val = nullptr;
	};
	void set(TTreeReader & ttreereader, std::string branchname, Bool_t _verbose = 1) {
		// if(!branchExists(branchname, ttreereader.GetTree())){
		// 	std::cout<<"\t\tError! Branch "<<branchname<<" does not exist in TTree "<<ttreereader.GetTree()->GetName()<<std::endl;
		// 	exit (EXIT_FAILURE);
		// }

		// if (!branchExists(branchname, ttreereader.GetTree())) {
		// 	std::cerr << "Error! Branch '" << branchname << "' not found!" << '\n';
		// 	exit(EXIT_FAILURE);
		// }

		val = new TTreeReaderArray<anytype>(ttreereader,branchname.c_str());

		if (_verbose) std::cout<<"\t\tInitialized array branch <"<< boost::typeindex::type_id<anytype>().pretty_name() <<"> "<<branchname<<std::endl;
	};
	anytype operator[](UInt_t index) {
		return (*(this->val)).At(index);
	};
	anytype at(UInt_t index) {
		return (*(this->val)).At(index);
	};
	UInt_t size() {
		return (*(this->val)).GetSize();
	};
	//returns pointer to 0th position
	operator anytype* () const {
		return &(*(*(this->val)).begin());
	};
};


template <typename anytype>
struct vector_association {
	std::vector<anytype> *vec = nullptr;
	anytype *var = nullptr;
	vector_association(anytype *_var, std::vector<anytype> *_vec): var(_var), vec(_vec) {};
	vector_association() {};
	~vector_association() {};
	anytype operator[](UInt_t index) {
		return vec->at(index);
	}
	void push_back() {
		vec->push_back(*var);
	};
	void clear() {
		vec->clear();
	};
};


struct plot_variable {
	Float_t *xptr = nullptr;
	Bool_t destroyPtr = 0;
	Float_t xmin;
	Float_t xmax;
	ULong_t nbins;
	Double_t *xBins = nullptr;
	std::string xtitle;
	std::string xunit;

	plot_variable(Float_t _xmin, Float_t _xmax, ULong_t _nbins,	std::string _xtitle = "", std::string _xunit = "") {
		set(_xmin, _xmax, _nbins, _xtitle, _xunit);
	};

	plot_variable(const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = "") {
		set(_xBins, _nbins, _xtitle, _xunit);
	};

	plot_variable(Float_t &_xptr, Float_t _xmin, Float_t _xmax, ULong_t _nbins,	std::string _xtitle = "", std::string _xunit = "") {
		set(_xptr, _xmin, _xmax, _nbins, _xtitle, _xunit);
	};

	plot_variable(Float_t &_xptr, const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = "") {
		set(_xptr, _xBins, _nbins, _xtitle, _xunit);
	};

	plot_variable() {};

	void set(Float_t _xmin, Float_t _xmax, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = "") {
		xptr = new Float_t;
		destroyPtr = 1;
		xmin = _xmin;
		xmax = _xmax;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
	};

	void set(const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = "") {
		xptr = new Float_t;
		destroyPtr = 1;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
		xBins = new Double_t[nbins+1];
		for (ULong_t i = 0; i <= _nbins; i++) {
			xBins[i] = _xBins[i];
		}
	};

	void set(Float_t &_xptr, Float_t _xmin, Float_t _xmax, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = "") {
		xptr = &_xptr;
		xmin = _xmin;
		xmax = _xmax;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
	};

	void set(Float_t &_xptr, const Double_t *_xBins, ULong_t _nbins, std::string _xtitle = "", std::string _xunit = "") {
		xptr = &_xptr;
		nbins = _nbins;
		xtitle = _xtitle;
		xunit = _xunit;
		xBins = new Double_t[nbins+1];
		for (ULong_t i = 0; i <= _nbins; i++) {
			xBins[i] = _xBins[i];
		}
	};

	~plot_variable() {
		if (destroyPtr) delete xptr;
		delete [] xBins;
	};

	operator Float_t () {
		return *xptr;
	};

	void operator = (const Float_t & assignVal) {
		*xptr = assignVal;
	};
};


struct histogram_template {
	const plot_variable *var = nullptr;
	const Float_t * varPtr = nullptr;
	TH1F* hist = nullptr;
	std::string histTitle;
	std::string histName;
	histogram_template(const plot_variable &_var, std::string _histTitle="", std::string _histName="", TH1F* _hist=nullptr) {
		set(_var, _histTitle, _histName, _hist);
	};

	histogram_template() {};

	void set(const plot_variable &_var, std::string _histTitle="", std::string _histName="", TH1F* _hist=nullptr) {
		var = &_var;
		varPtr = _var.xptr;
		histTitle = _histTitle;
		hist = _hist;
		if (hist) isUser = 1;
	};

	void initializehist(std::string name_prefix="", std::string title_prefix="", Bool_t Yunit = 1, TFile *_file = nullptr, std::string _TFileDir="") {
		if (isUser) std::cout<<"\tWarning: Data member .hist will not point to user-set histogram!"<<std::endl;
		std::string binwidth_string = "";
		if (var->xBins == nullptr && Yunit) {
			binwidth_string = removeTrailingZeros((var->xmax - var->xmin)/(Float_t)var->nbins);
			if (!binwidth_string.compare("1") && !var->xunit.empty()) binwidth_string.pop_back();
			else binwidth_string += " ";
			binwidth_string = "/"+binwidth_string;
			binwidth_string += var->xunit;
			rtrim(binwidth_string);
		}
		if (histName.empty()) histName = removeNonAlpha(name_prefix) + "_1D_" + removeNonAlpha(var->xtitle);
		trim(histName);
		trim(histTitle);
		std::string titles = title_prefix + histTitle + ";" + var->xtitle + (var->xunit.empty()?"":(" ["+var->xunit+"]")) + ";" + "# of events"+binwidth_string;
		if (var->xBins != nullptr) hist = new TH1F(histName.c_str(), titles.c_str(), var->nbins, var->xBins);
		else hist = new TH1F(histName.c_str(), titles.c_str(), var->nbins, var->xmin, var->xmax);
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->CenterTitle();
		hist->Sumw2();
		if (_file) {
			mkTFileDir(_file, _TFileDir);
			hist->SetDirectory(_file->GetDirectory(_TFileDir.c_str()));
		}
		isUser = 0;
		std::cout<<"\t\tInitialized TH1F "<<histName<<std::endl;
	};

	void fill(Double_t weight = 1.0) {
		if (!hist) {
			std::cout<<"Cannot fill! TH1F ("<<var->xtitle <<") is uninitialized!";
			return;
		}
		hist->Fill(*varPtr, weight);
	};

	~histogram_template() {
		// if(hist) delete hist;
	}
private:
	Bool_t isUser = 0;
};


struct twoDhistogram_template {
	const plot_variable *xvar = nullptr;
	const plot_variable *yvar = nullptr;

	const Float_t * xPtr = nullptr;
	const Float_t * yPtr = nullptr;

	TH2F* hist = nullptr;
	std::string histTitle;
	std::string histName;

	twoDhistogram_template(const plot_variable &_xvar, const plot_variable &_yvar, std::string _histTitle = "", std::string _histName="", TH2F* _hist = nullptr) {
		set(_xvar, _yvar, _histTitle, _histName, _hist);
	};

	twoDhistogram_template() {};

	void set(const plot_variable &_xvar, const plot_variable &_yvar, std::string _histTitle = "", std::string _histName="", TH2F* _hist = nullptr) {
		xvar = &_xvar;
		yvar = &_yvar;
		xPtr = _xvar.xptr;
		yPtr = _yvar.xptr;
		histTitle = _histTitle;
		hist = _hist;
		if (hist) isUser = 1;
	};

	void fill(Double_t weight = 1.0) {
		if (!hist) {
			std::cout<<"Cannot fill! TH2F ("<<xvar->xtitle << " VS "<< yvar->xtitle<<") is uninitialized!";
			return;
		}
		hist->Fill(*xPtr, *yPtr, weight);
	};

	void initializehist(std::string name_prefix = "", std::string title_prefix="", TFile *_file = nullptr, std::string _TFileDir="") {
		if (isUser) std::cout<<"\tWarning: Data member .hist will not point to user-set histogram!"<<std::endl;
		if (histName.empty())histName = removeNonAlpha(name_prefix) + "_2D_" + removeNonAlpha(xvar->xtitle +"\\ VS\\ " + yvar->xtitle);
		trim(histName);
		trim(histTitle);
		std::string titles = title_prefix + histTitle + ";" + xvar->xtitle + (xvar->xunit.empty()?"":(" ["+xvar->xunit+"]")) + ";" + yvar->xtitle + (yvar->xunit.empty()?"":(" ["+yvar->xunit+"]")) ;
		if (xvar->xBins != nullptr && yvar->xBins == nullptr) {
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xBins, yvar->nbins, yvar->xmin, yvar->xmax);
		} else if (xvar->xBins == nullptr && yvar->xBins != nullptr) {
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xmin, xvar->xmax, yvar->nbins, yvar->xBins);
		} else if (xvar->xBins != nullptr && yvar->xBins != nullptr) {
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xBins, yvar->nbins, yvar->xBins);
		} else {
			hist = new TH2F(histName.c_str(), titles.c_str(), xvar->nbins, xvar->xmin, xvar->xmax, yvar->nbins, yvar->xmin, yvar->xmax);
		}
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->CenterTitle();
		hist->SetTitle(histTitle.c_str());
		hist->Sumw2();

		if (_file) {
			mkTFileDir(_file, _TFileDir);
			hist->SetDirectory(_file->GetDirectory(_TFileDir.c_str()));
		}

		isUser = 0;
		std::cout<<"\t\tInitialized TH2F "<<histName<<std::endl;
	};

	~twoDhistogram_template() {
		// if(hist) delete hist;
	}
private:
	Bool_t isUser = 0;
};


struct bit_histogram_template: histogram_template {
	Bool_t *selector = nullptr;
	std::string prefix;
	bit_histogram_template() {};

	bit_histogram_template(const plot_variable &_var, Bool_t *_selector, std::string _prefix) {
		bitHistSet(_var, _selector, _prefix);
	};

	void bitHistSet(const plot_variable &_var, Bool_t *_selector, std::string _prefix) {
		selector = _selector;
		prefix = _prefix;
		set(_var);
	};

	void bitHistInit() {
		initializehist(prefix, prefix);
	};

	void fillBit(Double_t _weight=1.) {
		if (*selector) {
			fill(_weight);
		} else {
			hist->Fill(std::numeric_limits<Float_t>::max(), _weight);
		}
	};
};


struct bit_twoDhistogram_template: twoDhistogram_template {
	Bool_t *selector = nullptr;
	std::string prefix;
	bit_twoDhistogram_template() {};

	bit_twoDhistogram_template(const plot_variable &_xvar, const plot_variable &_yvar, Bool_t *_selector, std::string _prefix) {
		bitHistSet(_xvar, _yvar, _selector, _prefix);
	};

	void bitHistSet(const plot_variable &_xvar, const plot_variable &_yvar, Bool_t *_selector, std::string _prefix) {
		selector = _selector;
		prefix = _prefix;
		set(_xvar, _yvar);
	};

	void bitHistInit() {
		initializehist(prefix, prefix);
	};

	void fillBit(Double_t _weight=1.) {
		if (*selector) {
			fill(_weight);
		} else {
			hist->Fill(*(xvar->xptr), std::numeric_limits<Float_t>::max(), _weight);
		}
	};
};

struct signal_atts {
	signal_atts(std::string _couplingname, std::string _legend, std::string _color, Int_t _markerstyle) : couplingname(_couplingname), legend(_legend), color(_color), markerstyle(_markerstyle) {
	}

	signal_atts() {
	}
	std::string couplingname;
	std::string legend;
	std::string color;
	Int_t markerstyle;
	// std::string operator[]{
	// 	return couplingname;
	// };
};


struct sample {

	sample(std::string _ntuple, std::string _legend = "", Int_t _marker = 20, std::string _color = "#252525", Bool_t _drawLine = 0, TFile *_file = NULL, Float_t _luminosity=-999) {
		set(_ntuple, _legend, _marker, _color, _drawLine,_file, _luminosity);
	}

	sample() {
	}
	std::string ntuple;
	std::string legend;
	Int_t marker;
	Int_t lineStyle = 1;
	std::string color;
	Bool_t drawLine = 0;
	TFile *file;
	Float_t luminosity = 0.;

	void set(std::string _ntuple, std::string _legend = "", Int_t _marker = 20, std::string _color = "#252525", Bool_t _drawLine = 0, TFile *_file = NULL, Float_t _luminosity=-999) {
		ntuple= _ntuple;
		legend = _legend;
		marker = _marker;
		color = _color;
		drawLine = _drawLine;
		file = _file;
		luminosity = _luminosity;

		std::cout<<"\tInitialized sample: "<<std::endl<<
		"\t\tntuple: "<<ntuple<<std::endl<<
		"\t\tlegend: "<<legend<<std::endl<<
		"\t\tmarker: "<<marker<<std::endl<<
		"\t\tcolor: "<<color<<std::endl<<
		"\t\tdrawLine: "<<drawLine<<std::endl<<
		"\t\tfile: "<<file<<std::endl<<
		"\t\tluminosity: "<<luminosity<<std::endl;
	};

	void assignAtt(TH1 *_hist, Float_t _markerSize=1.5, Float_t _lineWidth=2.) {
		if (marker>0) {
			_hist->SetMarkerStyle(marker);
			_hist->SetMarkerSize(_markerSize);
			_hist->SetMarkerColor(TColor::GetColor(color.c_str()));
			_hist->SetLineColor(TColor::GetColor(color.c_str()));
			_hist->SetLineWidth(_lineWidth);
		} else {
			_hist->SetFillStyle((-marker));
			_hist->SetFillColor(TColor::GetColor(color.c_str()));
			_hist->SetLineColor(TColor::GetColor(color.c_str()));
			_hist->SetLineWidth(_lineWidth);
		}
		_hist->SetLineStyle(lineStyle);
	};

	void setLineStyle(Int_t _lineStyle) {
		lineStyle = _lineStyle;
	};
};


struct Profile2D {
	TH2D *hist = nullptr;
	std::vector<Double_t> _bin_entries;
	UInt_t nBinsTot;
	Profile2D(std::string _name, std::string _title, UInt_t _nBinsX, Double_t _xMin, Double_t _xMax, UInt_t _nBinsY, Double_t _yMin, Double_t _yMax) {
		set(_name, _title, _nBinsX, _xMin, _xMax, _nBinsY, _yMin, _yMax);
	};
	Profile2D() {};
	~Profile2D() {};
	void set(std::string _name, std::string _title, UInt_t _nBinsX, Double_t _xMin, Double_t _xMax, UInt_t _nBinsY, Double_t _yMin, Double_t _yMax) {
		hist = new TH2D(_name.c_str(), _title.c_str(), _nBinsX, _xMin, _xMax, _nBinsY, _yMin, _yMax);
		hist->GetXaxis()->CenterTitle();
		hist->GetYaxis()->CenterTitle();
		nBinsTot = (_nBinsX+2) * (_nBinsY+2);
		_bin_entries.reserve(nBinsTot);
		for (UInt_t i = 0; i < nBinsTot; i++) {
			_bin_entries.push_back(0.);
		}
	};
	void fill(Double_t _x, Double_t _y, Double_t _Zvalue, Double_t _weight=1.0) {
		hist->Fill(_x, _y, _Zvalue * _weight);
		UInt_t whichBin = hist->FindBin(_x, _y);
		_bin_entries[whichBin] += _weight;
	}
	TH2D * getProfile() {
		std::string profile_name = (std::string) hist->GetName() + "_profile";
		std::string profile_title = "Profile\\ " + (std::string) hist->GetTitle();
		TH2D *profile = (TH2D*) hist->Clone(profile_name.c_str());
		profile->Reset();
		profile->GetXaxis()->CenterTitle();
		profile->GetYaxis()->CenterTitle();
		profile->SetTitle(profile_title.c_str());
		Double_t sumEntries = std::accumulate(_bin_entries.begin(), _bin_entries.end(), 0);
		for (UInt_t i = 0; i < nBinsTot; i++) {
			Double_t _mean = _bin_entries[i] > 0. ? hist->GetBinContent(i)/_bin_entries[i] : 0. ;
			if (_bin_entries[i] > 0.) profile->SetBinContent(i, _mean);
		}
		return profile;
	}
	operator TH2D* () {
		return hist;
	}
};

template <typename anytype>
struct indexer {
	indexer() {};
	indexer(std::vector<anytype> _list) {
		init(_list);
	};
	~indexer() {};

	void init(std::vector<anytype> _list) {
		for (UInt_t iIndex = 0; iIndex < _list.size(); iIndex++) {
			theIndex[_list[iIndex]]	= iIndex;
		}
	};

	UInt_t size() {
		return theIndex.size();
	};

	std::map<anytype, Int_t> theIndex;

	Int_t operator[](anytype needle) {
		if ( theIndex.find(needle) != theIndex.end()) return theIndex[needle];
		else return -999;
	};
};

class CSVReader {
private:

	void parseFile(std::string fileName, std::string delimeter = ",") {
		if (!file_exists(fileName)) {
			std::cout<<"Error! Cannot CSV-parse! File does not exist:"<<fileName<<std::endl;
			return;
		}
		dataList.clear();
		std::ifstream file(fileName);
		std::string line = "";
		while (getline(file, line)) {
			trim(line);
			std::vector<std::string> line_split = split_string(line, delimeter);
			if (line_split.empty()) continue;
			dataList.push_back(line_split);
		}
		file.close();
	};

public:
	CSVReader() {};

	CSVReader(std::string fileName, std::string delimeter = ",") {
		parseFile(fileName, delimeter);
	};

	std::vector<std::vector<std::string>> getData() {
		return dataList;
	};

	std::vector<std::vector<std::string>> dataList;

	std::vector<std::string> operator[](UInt_t index) {
		if (dataList.empty()) return std::vector<std::string> {};
		return dataList[index];
	};

	operator std::vector<std::vector<std::string>> () const {
		return dataList;
	};

	UInt_t size() {
		return dataList.size();
	};
};


struct kFactorMap {
	std::map<Float_t, Float_t> theMap;

	Float_t xMax;
	Float_t xMin;

	Bool_t inValid = 1;

	kFactorMap() {};

	kFactorMap(std::string histFile, std::string histName, Bool_t verbose = 1) {
		init(histFile, histName, verbose);
	};

	void init(std::string histFile, std::string histName, Bool_t verbose = 1) {

		TH1* 	KFactorHist = getHistFromFile(histName, histFile, verbose);

		KFactorHist->SetName(randStr().c_str());

		if (!KFactorHist) {
			std::cout<<"Warning! Kfactor histogram "<<histName<<" not found in file "<<histFile<<
			"! Will return 1 for all queries."<<std::endl;
			return;
		}

		KFactorHist->SetName(randStr(6).c_str());

		xMin = KFactorHist->GetXaxis()->GetBinCenter(1);
		xMax = KFactorHist->GetXaxis()->GetBinCenter(KFactorHist->GetNbinsX());

		for (Int_t iBin = 1; iBin <= KFactorHist->GetNbinsX(); iBin++) {



			Float_t iBinMax 		=		KFactorHist->GetXaxis()->GetBinUpEdge(iBin);
			Float_t iBinKfac 		=		KFactorHist->GetBinContent(iBin);
			theMap[iBinMax]			=		iBinKfac;

			if (verbose) {
				Float_t iBinCenter 		=		KFactorHist->GetXaxis()->GetBinCenter(iBin);
				std::cout<<KFactorHist->GetXaxis()->GetBinLowEdge(iBin)<<"-"<<iBinMax<<"\t"<<iBinKfac<<std::endl;
			}
		}

		delete KFactorHist;

		if (theMap.empty()) {
			inValid = 1;
			std::cout<<"Warning! Kfactor map taken from "<<histName<<" in file "<<histFile<<
			"is empty!"<<std::endl;
		}

		if (verbose) std::cout<<"kfactors loaded from TH1 "<<histName<< " in "<<histFile <<std::endl;
	};

	Float_t getKFac(const Float_t pTval) {

		if (theMap.empty()) return 1.;

		std::map<Float_t, Float_t>::iterator iBin;

		Float_t searchPfVal = pTval;
		searchPfVal = std::min(xMax, searchPfVal);
		searchPfVal = std::max(xMin, searchPfVal);

		iBin = theMap.lower_bound(searchPfVal);

		return iBin->second;
	};

};


struct isoCorrMap {

	std::map<Float_t, TF1*> theMap;
	std::map<Float_t, Float_t> zeroVals;
	std::string mapFile;

	isoCorrMap() {};

	isoCorrMap(std::string mapFile, UInt_t fColumn = 3, Bool_t _verbose=1, std::string _delimiter=";", Int_t _firstLine = 0) {
		init(mapFile, fColumn, _verbose, _delimiter, _firstLine);
	};

	~isoCorrMap() {
		for (std::map<Float_t, TF1*>::iterator iEl 		= 	theMap.begin(); iEl != theMap.end(); iEl++) {
			delete iEl->second;
			// iEl->second = nullptr;
		}

		theMap.clear();
	};

	void init(std::string _mapFile, UInt_t fColumn = 3, Bool_t _verbose=1, std::string _delimiter=";", Int_t _firstLine = 0) {
		if (!file_exists(_mapFile)) {
			std::cout<<"\tError! Isolation map file "<< _mapFile<<" does not exist! "<<mapFile<<std::endl;
			return;
		}

		mapFile = _mapFile;
		CSVReader readFile(mapFile, _delimiter);
		std::vector<std::vector<std::string>> isCorrData 		= 	readFile.getData();

		std::cout<<"Loading isolation corrections from "<<_mapFile<<std::endl;

		for (UInt_t iRow = _firstLine; iRow < isCorrData.size(); iRow++) {

			if (isCorrData[iRow].size() < fColumn + 1) continue;

			Float_t uBound 										= 	std::stof(isCorrData[iRow][1]);

			if (_verbose) {
				std::cout<<"\t\t\t"<<uBound<<"\t\t"<<isCorrData[iRow][fColumn]<<std::endl;
			}

			TF1* isCorrFunc 									= 	new TF1(removeNonAlpha(getFileName(mapFile) + isCorrData[iRow][1]).c_str(), isCorrData[iRow][fColumn].c_str(), 0.,
				13000.);
			theMap[uBound] 										= 	isCorrFunc;

			zeroVals[uBound]									=	isCorrFunc->Eval(0);
		}

	};

	Float_t getIsoCorr(const Float_t absEta, const Float_t sendaryEn, Bool_t _keepZero = 0) {
		if (theMap.empty()) return 0.;
		std::map<Float_t, TF1*>::iterator iBin = theMap.lower_bound(absEta);
		std::map<Float_t, Float_t>::iterator iBinZero = zeroVals.lower_bound(absEta);

		if (iBin == theMap.end()) {
			std::cout<<"Error! Eta "<<absEta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		if (_keepZero) return iBin->second->Eval(sendaryEn);
		return (iBin->second->Eval(sendaryEn) - iBinZero->second);
	};

	Float_t getEffectiveAreaAbs(const Float_t eta, const Float_t sendaryEn, Bool_t _keepZero = 0) {
		if (theMap.empty()) return 0.;
		Float_t 		absEta = std::abs(eta);
		std::map<Float_t, TF1*>::iterator iBin = theMap.lower_bound(absEta);
		std::map<Float_t, Float_t>::iterator iBinZero = zeroVals.lower_bound(absEta);

		if (iBin == theMap.end()) {
			std::cout<<"Error! Eta "<<eta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		if (_keepZero) return iBin->second->Eval(sendaryEn);
		return (iBin->second->Eval(sendaryEn) - iBinZero->second);

	};
};


struct effectiveAreaMap {
	std::map<Float_t, Float_t> effectiveAreas;
	std::string mapFile;

	effectiveAreaMap() {};

	effectiveAreaMap(std::string mapFile, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2) {
		init(mapFile, _verbose, _delimiter, _firstLine);
	};

	void init(std::string _mapFile, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2) {
		if (!file_exists(_mapFile)) {
			std::cout<<"\tError! File does not exist! "<<mapFile<<std::endl;
			exit(EXIT_FAILURE);
		}

		mapFile = _mapFile;
		CSVReader readFile(mapFile, _delimiter);
		std::vector<std::vector<std::string>> effAreaData 		= 	readFile.getData();

		for (UInt_t iRow = _firstLine; iRow < effAreaData.size(); iRow++) {
			Float_t uBound 										= 	std::stof(effAreaData[iRow][1]);
			Float_t effArea 									= 	std::stof(effAreaData[iRow][2]);
			effectiveAreas[uBound] 								= 	effArea;
		}

		if (_verbose) {
			std::cout<<"\tEffective areas loaded from file: "<< mapFile<<std::endl;
			std::cout<<"\t\t\t|eta|\t\tEA:"<<std::endl;
			for (std::map<Float_t, Float_t>::iterator iEl 		= 	effectiveAreas.begin(); iEl != effectiveAreas.end(); iEl++) {
				std::cout<<"\t\t\t"<<iEl->first<<"\t\t"<<iEl->second<<std::endl;
			}
		}
	};

	Float_t getEffectiveArea(Float_t absEta) {
		if (effectiveAreas.empty()) return 0.;
		std::map<Float_t, Float_t>::iterator iBin = effectiveAreas.lower_bound(absEta);
		if (iBin == effectiveAreas.end()) {
			// std::cout<<"Error! Eta "<<absEta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		return iBin->second;
	};

	Float_t getEffectiveAreaAbs(Float_t eta) {
		if (effectiveAreas.empty()) return 0.;
		std::map<Float_t, Float_t>::iterator iBin = effectiveAreas.lower_bound(std::abs(eta));
		if (iBin == effectiveAreas.end()) {
			// std::cout<<"Error! Eta "<<eta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		return iBin->second;
	};
};

struct isoPtScalingMap {

	std::map<Float_t, std::pair<Float_t,Float_t>> ptScalingCoeffs;
	std::map<Float_t, std::pair<Float_t,Float_t>>::iterator etaBin;
	Bool_t isQuadratic;
	std::string mapFile;

	isoPtScalingMap() {};

	isoPtScalingMap(std::string _mapFile, Bool_t _isQuadratic, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2) {
		init(_mapFile, _isQuadratic, _verbose, _delimiter, _firstLine);
	};

	void init(std::string _mapFile, Bool_t _isQuadratic, Bool_t _verbose=1, std::string _delimiter="        ", Int_t _firstLine = 2) {
		if (!file_exists(_mapFile)) {
			std::cout<<"\tError! Isolation scaling file does not exist! "<<_mapFile<<std::endl;
			exit(EXIT_FAILURE);
		}
		mapFile = _mapFile;
		isQuadratic = _isQuadratic;

		CSVReader readFile(_mapFile, _delimiter);
		std::vector<std::vector<std::string>> ptScalingData 		= 	readFile.getData();

		for (UInt_t iRow = _firstLine; iRow < ptScalingData.size(); iRow++) {
			Float_t uBound 										= 	std::stof(ptScalingData[iRow][1]);
			Float_t c1 											= 	std::stof(ptScalingData[iRow][3]);
			Float_t c2 											= 	_isQuadratic ? std::stof(ptScalingData[iRow][4]) : 0.;
			ptScalingCoeffs[uBound] 								= 	std::make_pair(c1,c2);
		}

		if (_verbose) {
			std::cout<<"\tIsolation pT scalings loaded from file: "<< _mapFile<<std::endl;
			std::cout<<"\t\t\t|eta|\t\t\tc1:\t\t\tc2"<<std::endl;
			for (std::map<Float_t, std::pair<Float_t,Float_t>>::iterator iEl	= 	ptScalingCoeffs.begin(); iEl != ptScalingCoeffs.end(); iEl++) {
				std::cout<<"\t\t\t"<<iEl->first<<"\t\t\t"<<iEl->second.first<<"\t\t"<<iEl->second.second<<std::endl;
			}
		}
	};

	Float_t getPtScaling(Float_t _absEta, Float_t _pT) {
		etaBin = ptScalingCoeffs.lower_bound(_absEta);

		if (etaBin == ptScalingCoeffs.end()) {
			// std::cout<<"Error! Eta "<< _absEta <<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}

		Float_t scalingCorrection = (etaBin->second.first)*_pT + (isQuadratic ? (etaBin->second.second)*_pT*_pT : 0.);
		return scalingCorrection;
	};

	Float_t getPtScalingAbs(Float_t _eta, Float_t _pT) {
		etaBin = ptScalingCoeffs.lower_bound(std::abs(_eta));

		if (etaBin == ptScalingCoeffs.end()) {
			std::cout<<"Error! Eta "<< _eta <<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}

		Float_t scalingCorrection = (etaBin->second.first)*_pT + (isQuadratic ? (etaBin->second.second)*_pT*_pT : 0.);
		return scalingCorrection;
	};

};


struct coronaCorrections {
	std::vector<std::vector<Float_t>>							percentiles;
	std::vector<std::vector<std::vector<Float_t>>> 				percentileLines;
	std::vector<Float_t> 										etaLimits;


	coronaCorrections(std::string mapFile, std::string pDelim = ";", std::string sDelim=",", Bool_t verbose = 1) {
		init(mapFile, pDelim, sDelim, verbose);
	};

	coronaCorrections() {};

	void init(std::string mapFile, std::string pDelim = ";", std::string sDelim=",", Bool_t verbose = 1) {
		percentileLines.clear();
		etaLimits.clear();
		percentiles.clear();

		CSVReader readFile(mapFile, pDelim);
		std::vector<std::vector<std::string>> effAreaData 		= 	readFile.getData();

		if (verbose) std::cout<<"\tLoading Corona corrections from file "<<mapFile<<std::endl;

		for (UInt_t i = 0; i < effAreaData.size(); i++) {

			etaLimits.push_back(strToFloatList(effAreaData[i][0], sDelim)[1]);

			std::vector<std::vector<Float_t>>	etaLines;
			std::vector<Float_t>				iPercentiles;

			for (UInt_t j = 1; j < effAreaData[i].size(); j++) {
				std::vector<Float_t> lineDef =	strToFloatList(effAreaData[i][j], sDelim);

				if (lineDef.size() != 3) continue;

				iPercentiles.push_back(lineDef[0]);
				etaLines.push_back({lineDef[1], lineDef[2]});
			};

			percentileLines.push_back(etaLines);
			percentiles.push_back(iPercentiles);
		};

		if (verbose) printData();
	};

	void printData() {
		for (UInt_t iEta = 0; iEta< etaLimits.size(); iEta++) {
			std::cout<<"eta < "<<etaLimits[iEta]<<std::endl;
			for (UInt_t iP = 0; iP < percentileLines[iEta].size(); iP++) {
				std::cout<<"\t"<<percentiles[iEta][iP]*100.<<"%-tile"<<std::endl;
				for (UInt_t iData = 0; iData < percentileLines[iEta][iP].size(); iData++) {
					std::cout<<"\t\t"<<percentileLines[iEta][iP][iData]<<"\t";
				}
				std::cout<<std::endl;
			}
		}
	};

	Float_t getCorrection(Float_t _absEta, Float_t _iso, Float_t _rho, Bool_t verbose = 0) {
		Int_t etaIndex = -999;
		for (UInt_t iEta = 0; iEta < etaLimits.size(); iEta++) {
			if (_absEta < etaLimits[iEta]) {
				etaIndex = iEta;
				break;
			};
		};

		if (etaIndex < 0) return 0.;

		Bool_t 	foundCorrection 		= 		0;
		Int_t 	corrP 					= 		0;

		for (Int_t iPercentile = percentileLines[etaIndex].size() - 1; iPercentile > -1 ; iPercentile--) {
			Float_t yLine 				=		percentileLines[etaIndex][iPercentile][0]*_rho + percentileLines[etaIndex][iPercentile][1];
			if (verbose) std::cout<<percentiles[etaIndex][iPercentile]*100<<"%-tile \t yLine = "<<yLine<<std::endl;
			if (_iso < yLine) {
				corrP 			= 		iPercentile;
				break;
			};
		}

		if (verbose) std::cout<<"eta = "<< etaLimits[etaIndex]<<"\t % = "<< percentiles[etaIndex][corrP]*100<<"\t m = "<<percentileLines[etaIndex][corrP][0]<<"\t c = "<<percentileLines[etaIndex][corrP][1]<<std::endl;

		Float_t corrrection 				= 		percentileLines[etaIndex][corrP][0]*_rho;

		return corrrection;
	};

};


struct parseOptions {

	std::string delimeter;
	std::string optFile;
	std::map<std::string, std::string> optMap;
	Bool_t isParsed = false;

	// static std::string addNspaces(std::string _inStr, UInt_t _nSpace, Bool_t _pre = 0) {
	// 	std::string spaceStrToAdd="";
	// 	for (UInt_t i = 0; i <_nSpace; i++) {
	// 		spaceStrToAdd += " ";
	// 	}
	// 	std::string resultStr = _pre ? (spaceStrToAdd + _inStr) : (_inStr + spaceStrToAdd);
	// 	return resultStr;
	// }

	void writeCard(std::string _writePath) {
		std::string writeStr =  "#### parseOptions inputs from\n#### " + optFile +"\n#### @ " + getCurrentTime() +"\n";
		for (std::pair<std::string, std::string> it : optMap) {
			writeStr += it.first + " \t\t\t" + delimeter + "\t\t\t " + it.second + "\n\n";
		}
		if (_writePath.empty()) std::cout<<writeStr<<std::endl;
		else writeStringToFile(writeStr, _writePath, 0);
	}

	// void writeCard(std::string _writePath, std::string _auxSeparator = ";") {
	// 	if (_auxSeparator == delimeter) {
	// 		std::cout<<__LINE__<<" parseOptions error! Auxiliary separator ["<<_auxSeparator<<"] is identical to primary ["<<delimeter<<"]"<<std::endl;
	// 		return;
	// 	}
	// 	UInt_t maxOptWidth = 20;
	// 	for (std::pair<std::string, std::string> it : optMap) {
	// 		if (it.first.length() > maxOptWidth) maxOptWidth = it.first.length();
	// 	}
	// 	std::string writeStr =  "#### parseOptions inputs from\n#### " + optFile +"\n#### @ " + getCurrentTime() +"\n";
	// 	for (UInt_t i = 0; i < maxOptWidth; i++) {
	// 		writeStr +="#";
	// 	}
	// 	writeStr +="\n";
	// 	UInt_t lenAfterPrimDelim = maxOptWidth + delimeter.length() + 2;
	// 	for (std::pair<std::string, std::string> it : optMap) {
	// 		// std::string keyStr = (it.first.length() < maxOptWidth) ? addNspaces(it.first, maxOptWidth - it.first.length()) : it.first;
	// 		// writeStr += keyStr + " " + delimeter + " ";
	// 		writeStr += it.first + " " + delimeter + " ";
	// 		std::vector<std::string> keySubOpts = split_string(it.second, _auxSeparator);
	// 		if (keySubOpts.size() == 1) writeStr += keySubOpts[0] + "\n";
	// 		else continue;
	// 		for (std::string iKeySubOpt : keySubOpts) {
	// 			// writeStr = addNspaces(writeStr, lenAfterPrimDelim, 0);
	// 			writeStr += iKeySubOpt + _auxSeparator + "\n";
	// 		}
	// 		UInt_t nStrToRm = _auxSeparator.length() + 2;
	// 		for (UInt_t i = 0; i < nStrToRm; i++) {
	// 			writeStr.pop_back();
	// 		}
	// 		writeStr += "\n";
	// 	}

	// 	if (_writePath.empty()) std::cout<<writeStr<<std::endl;
	// 	else writeStringToFile(writeStr, _writePath, 0);
	// }

	parseOptions(std::string _optFile, std::string _delimiter=",", Bool_t _verbose = 0) {
		optFile = _optFile;
		parseIt(optFile, _delimiter, _verbose);

	};

	parseOptions() {};

	Bool_t keyExists(std::string _key) {
		if (optMap.find(_key) == optMap.end()) {
			return 0;
		} else {
			return 1;
		}
	}

	void parseIt(std::string _optFile, std::string _delimiter=",", Bool_t _verbose = 0) {
		delimeter = _delimiter;
		optFile = _optFile;
		if (!file_exists(_optFile)) {
			std::cout<<"Error! Options file "<<_optFile<<" not found!"<<std::endl;
		}
		optMap.clear();
		if (_verbose) std::cout<<"\t\tOptions parsed from file "<<_optFile<<"... "<<std::endl<<std::endl;
		CSVReader _csvFile(_optFile, _delimiter);
		std::vector<std::vector<std::string>> _data = _csvFile.getData();

		std::string lastOption="";
		for (UInt_t i = 0; i < _data.size(); i++) {
			std::string _optName = _data[i][0];
			trim(_optName);
			if (_optName.empty()) continue;
			if (match("#*", _optName)) continue;
			if ((_data[i].size()==1)) {
				optMap[lastOption] += trim_copy(_data[i][0]);
				continue;
			};
			std::string _optVal = trim_copy(_data[i][1]);
			optMap[_optName] = _optVal;
			lastOption = _optName;
		}

		isParsed = true;

		if (_verbose) {
			for (std::map<std::string, std::string>::iterator optKey = optMap.begin(); optKey != optMap.end(); optKey++) {
				std::cout<<"\t\t\t"<<optKey->first<<"\t\t"<<optKey->second<<std::endl;
			}

			std::cout<<"\tOptions parsed @ "<<getCurrentTime()<<std::endl;
		};
	};

	Float_t getFloat(std::string _opt) {
		return std::stof(get(_opt));
	};

	Double_t getDouble(std::string _opt) {
		return std::stod(get(_opt));
	};

	Int_t getInt(std::string _opt) {
		return std::stoi(get(_opt));
	};

	Bool_t getStrictBool(std::string _opt) {
		Int_t optDecision = std::stoi(get(_opt));
		if (optDecision == 1) return true;
		return false;
	};

	std::map<std::string, std::string> getMap(std::string _opt, std::string _delimiter1, std::string _delimiter2) {

		std::vector<std::string> optMapElements = getList(_opt, _delimiter1);
		std::map<std::string, std::string> optSubMap;
		for (std::string iElement : optMapElements) {
			std::vector<std::string> elementSplit = split_string(iElement, _delimiter2);
			if (elementSplit.size() < 2) {
				std::cout<<"Error!  parseOptions::getMap() failed for key "<<_opt<<" with delimiters ["<<_delimiter1<<"] and ["<<_delimiter2<<"]"<<std::endl;
				return optSubMap;
			}
			optSubMap[elementSplit[0]] = elementSplit[1];
		}
		return optSubMap;
	}

// 	std::string getSubOpt(std::string _opt0, std::string _opt1, std::string _subDelimiter0, std::string _subDelimiter1) {
// 		if (!keyExists(_opt)) {
// 			cout<<"Error! Key \""<<_opt<<"\" not found in file "<<optFile <<" !"<<std::endl;
// 			return std::string("");
// 		}
// 		std::vector<std::string> subOpts = split_string(optMap.at(_opt), _subDelimiter0);
// 		std::string subOptVal = "";
// 		for(std::string iSubOpt : subOpts){
// 			std::vector<std::string> iSubOpt = split_string(iSubOpt, _subDelimiter1);
// 			if(iSubOpt[0] == _opt1){
// subOptVal =
// 			}
// 		}
// 		return optMap.at(_opt);
// 	};

	std::string get(std::string _opt) {
		if (!keyExists(_opt)) {
			cout<<"Key \""<<_opt<<"\" not found in file "<<optFile <<" !"<<std::endl;
			return std::string("");
		}
		return optMap.at(_opt);
	};

	void addOpt(std::string _OptKey, std::string _OptVal) {
		if (_OptKey.empty()) {
			std::cout<<"Warning: empty option key provided!"<<std::endl;
		} else {
			optMap[_OptKey] = _OptVal;
		}
	};

	const char* getCSTR(std::string _opt) {
		if (!keyExists(_opt)) {
			cout<<"Error! Key "<<_opt<<" not found!"<<std::endl;
		}
		return optMap.at(_opt).c_str();
	};

	std::string getListElement(std::string _opt, UShort_t iElement, std::string _delimiter = ",") {
		std::vector<std::string> theList = getList(_opt, _delimiter);
		if (theList.size() <= iElement) {
			std::cout<<"parseOptions error (file="<<optFile<<") :\n\t\tgetListElement(opt=" <<_opt<<", iElement ="<<iElement<<", delimiter = "<<_delimiter <<") called\n\t\twhile list size is "<<theList.size()<<std::endl;
			return "";
		}
		return theList[iElement];
	};

	std::vector<std::string> getList(std::string _opt, std::string _delimiter = ",") {
		std::string listStr = get(_opt);
		if (listStr.empty()) return std::vector<std::string>({});

		std::vector<std::string> theList = split_string(listStr, _delimiter, 1);
		if (theList.back().empty()) theList.pop_back();
		return theList;
	};

	std::vector<Float_t> getFloatList(std::string _opt, std::string _delimiter = ",") {
		return  strToFloatList(get(_opt), _delimiter);
	};

	std::vector<Double_t> getDoubleList(std::string _opt, std::string _delimiter = ",") {
		return  strToDoubleList(get(_opt), _delimiter);
	};

	std::vector<Int_t> getIntList(std::string _opt, std::string _delimiter = ",") {
		return  strToIntList(get(_opt), _delimiter);
	};

	std::vector<Int_t> getTColListFromHexList(std::string _opt, std::string _delimiter = ",") {

		std::vector<std::string> 	hexList 		=	split_string(get(_opt), _delimiter);
		std::vector<Int_t> 			tColList;

		for (std::string & hexCol : hexList) {
			tColList.push_back(hex2rootColor(hexCol));
		};

		return  tColList;
	};

	Int_t getTColFromHex(std::string _opt) {
		return hex2rootColor(get(_opt));
	};

	std::string & operator[](std::string _opt) {
		return optMap[_opt];
	};
};


struct bitBar2D {
	TH2F hist;
	std::vector<const Bool_t*> xVars;
	const Float_t *yVar = nullptr;

	bitBar2D() {};
	bitBar2D(std::vector<std::pair<Bool_t*, std::string>> _xVars, const plot_variable & _yVar) {
		init(_xVars, _yVar);
	};

	void init(std::vector<std::pair<Bool_t*, std::string>> _xVars, const plot_variable & _yVar) {
		yVar = _yVar.xptr;
		Int_t nxBins = _xVars.size();

		std::string histName = removeNonAlpha(_yVar.xtitle) + "__vs__";
		std::string histTitle = _yVar.xtitle + " vs (";
		for (UInt_t i = 0; i < _xVars.size(); i++) {
			histName += removeNonAlpha(_xVars[i].second)+"_";
			histTitle += _xVars[i].second + ", ";
			xVars.push_back(_xVars[i].first);
		}
		histName.pop_back();
		histTitle.pop_back();
		histTitle.pop_back();
		histTitle += ")";

		histTitle += ";;" + _yVar.xtitle;
		histTitle += _yVar.xunit.empty() ? "" : ("[" + _yVar.xunit + "]");
		if (_yVar.xBins == nullptr) hist = TH2F(histName.c_str(), histTitle.c_str(), nxBins, 0., (Float_t)nxBins, _yVar.nbins, _yVar.xmin, _yVar.xmax);
		else hist = TH2F(histName.c_str(), histTitle.c_str(), nxBins, 0., (Float_t)nxBins, _yVar.nbins, _yVar.xBins);

		for (UInt_t i = 0; i < _xVars.size(); i++) {
			hist.GetXaxis()->SetBinLabel(i+1, _xVars[i].second.c_str());
		}

		hist.GetYaxis()->CenterTitle();

		std::cout<<"\t\t\tInitialized bitBar2D "<<histTitle<<std::endl;
	}

	void fill(Float_t weight = 1.) {
		for (UInt_t i = 0; i < xVars.size(); i++) {
			hist.Fill((Float_t)*xVars[i]+0.01, *yVar, (Float_t)*xVars[i] * weight);
		}
	}
};

struct BGset {
	BGset() {};

	BGset(std::string _filenames, std::string _legend, std::string _hexColor, std::string _separator="", std::string _xsecFacator="1") {
		init(_filenames, _legend, _hexColor, _separator, _xsecFacator);
	};

	std::vector<std::string> fileNames;
	std::string legend;
	std::string hexColor;
	Double_t xsecFacator;

	void init(std::string _filenames, std::string _legend, std::string _hexColor, std::string _separator=",", std::string _xsecFacator="1") {
		fileNames = split_string(_filenames, _separator);
		legend = _legend;
		hexColor = _hexColor;
		trim(_xsecFacator);
		xsecFacator=std::stod(_xsecFacator);

		std::cout<<"\tInitialized background set:"<<std::endl<<"\t\tFiles:"<<std::endl;;

		for (std::string fileName : fileNames) {
			std::cout<<"\t\t\t"<<fileName<<std::endl;
		}

		std::cout<<"\t\tLegend: "<<legend<<std::endl;
		std::cout<<"\t\tColor: "<<hexColor<<std::endl;
		std::cout<<"\t\txSec k-factor: "<<xsecFacator<<std::endl;
	};

	void assignAttFill(TH1 *_hist, TLegend *_legend = nullptr) {
		_hist->SetFillColor(TColor::GetColor(hexColor.c_str()));
		if (_legend) {
			TLegendEntry *legEntry = _legend->AddEntry(_hist, legend.c_str(), "F");
			legEntry->SetTextColor(TColor::GetColor(hexColor.c_str()));
		}
	};

	void assignAttLine(TH1 *_hist, TLegend *_legend = nullptr) {
		_hist->SetLineColor(TColor::GetColor(hexColor.c_str()));
		if (_legend) {
			TLegendEntry *legEntry = _legend->AddEntry(_hist, legend.c_str(), "LFPE");
			legEntry->SetTextColor(TColor::GetColor(hexColor.c_str()));
		}
	};
};


template <typename anytype>
struct smarter_ptr {
	anytype * thePtr = nullptr;

	smarter_ptr() {};

	smarter_ptr(anytype *_thePtr):thePtr(_thePtr) {};

	~smarter_ptr() {
		delete thePtr;
		thePtr = nullptr;
	};

	anytype * get() {
		return thePtr;
	};

	anytype* operator->() {
		return thePtr;
	}

	smarter_ptr & operator = (const anytype * & _assign) {
		thePtr = _assign;
	};

	operator anytype*() {
		return thePtr;
	};

	// anytype operator -> anytype(){
	// 	return *thePtr;
	// };
};


struct scaleFactor1D {
//// https://indico.cern.ch/event/755652/contributions/3131649/attachments/1712933/2762198/Electron_ID_V2_SFs_on_2017_and_2018_data_20180910.pdf#search=Dalmin%20Pai
//// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
	scaleFactor1D() {};

	scaleFactor1D(std::string _histFilePath, std::string _histName) {
		init(_histFilePath, _histName);
	};

	scaleFactor1D(std::string _histFilePathHistNameStr) {
		init(_histFilePathHistNameStr);
	};

	~scaleFactor1D() {
		delete sfMap;
	};

	TH1* sfMap = nullptr;

	Float_t xMin;
	Float_t xMax;
	Float_t defVal;


	void init(std::string _histFilePath, std::string _histName, Float_t _defVal = 1.) {
		defVal = _defVal;

		sfMap = (TH1*) getHistFromFile(_histName, _histFilePath);

		if (!sfMap) {
			std::cout<<"Warning: "<<_histName<<" could not be loaded from "<<_histFilePath<<"! Will return "<<defVal<<" for SF!"<<std::endl;
			return;
		}

		std::cout<<"Loading scale factors..."<<std::endl;

		xMin = sfMap->GetXaxis()->GetBinCenter(1);
		xMax = sfMap->GetXaxis()->GetBinCenter(sfMap->GetNbinsX());
	};

	void init(std::string _histFilePathHistNameStr) {
		init(split_string(_histFilePathHistNameStr, ",")[0], split_string(_histFilePathHistNameStr, ",")[1]);
	};

	Float_t getSF(Float_t _var) {

		if (!sfMap) {
			return defVal;
		}

		_var = std::max(_var, xMin);
		_var = std::min(_var, xMax);

		Int_t bin = sfMap->FindBin(_var);
		return sfMap->GetBinContent(bin);
	};
};

struct pTEtaScaleFactor2D {
//// https://indico.cern.ch/event/755652/contributions/3131649/attachments/1712933/2762198/Electron_ID_V2_SFs_on_2017_and_2018_data_20180910.pdf#search=Dalmin%20Pai
//// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
	pTEtaScaleFactor2D() {};

	pTEtaScaleFactor2D(std::string _histFilePath, std::string _histName) {
		init(_histFilePath, _histName);
	};

	pTEtaScaleFactor2D(std::string _histFilePathHistNameStr) {
		init(_histFilePathHistNameStr);
	};

	~pTEtaScaleFactor2D() {
		delete sfMap;
	};

	TH2* sfMap = nullptr;
	Float_t ptMin;
	Float_t ptMax;


	void init(std::string _histFilePath, std::string _histName) {

		sfMap = (TH2*) getHistFromFile(_histName, _histFilePath);

		if (!sfMap) {
			std::cout<<"Warning: "<<_histName<<" could not be loaded from "<<_histFilePath<<"! Will return 1 for SF!"<<std::endl;
			return;
		}
		std::cout<<"Loading scale factors..."<<std::endl;
		ptMin = sfMap->GetXaxis()->GetBinCenter(1);
		ptMax = sfMap->GetXaxis()->GetBinCenter(sfMap->GetNbinsX());
	};

	void init(std::string _histFilePathHistNameStr) {
		init(split_string(_histFilePathHistNameStr, ",")[0], split_string(_histFilePathHistNameStr, ",")[1]);
	};

	Float_t getSF(Float_t _pT, Float_t _eta) {

		// if (!sfMap) return 1;

		_pT = std::max(_pT, ptMin);
		_pT = std::min(_pT, ptMax);

		Int_t bin = sfMap->FindBin(_pT, _eta);
		return sfMap->GetBinContent(bin);
	};
};


struct etaPtScaleFactor2D {

	//// https://indico.cern.ch/event/755652/contributions/3131649/attachments/1712933/2762198/Electron_ID_V2_SFs_on_2017_and_2018_data_20180910.pdf#search=Dalmin%20Pai
	//// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017

	etaPtScaleFactor2D() {};

	etaPtScaleFactor2D(std::string _histFilePath, std::string _histName) {
		init(_histFilePath, _histName);
	};

	etaPtScaleFactor2D(std::string _histFilePathHistNameStr) {
		init(_histFilePathHistNameStr);
	};

	~etaPtScaleFactor2D() {
		delete sfMap;
	};

	TH2* sfMap = nullptr;
	Float_t ptMin;
	Float_t ptMax;


	void init(std::string _histFilePath, std::string _histName) {

		sfMap = (TH2*) getHistFromFile(_histName, _histFilePath);
		if (!sfMap) {
			std::cout<<"Warning: "<<_histName<<" could not be loaded from "<<_histFilePath<<"! Will return 1 for SF!"<<std::endl;
			return;
		}

		std::cout<<"Loading scale factors..."<<std::endl;
		ptMin = sfMap->GetYaxis()->GetBinCenter(1);
		ptMax = sfMap->GetYaxis()->GetBinCenter(sfMap->GetNbinsY());
	};

	void init(std::string _histFilePathHistNameStr) {
		init(split_string(_histFilePathHistNameStr, ",")[0], split_string(_histFilePathHistNameStr, ",")[1]);
	};

	Float_t getSF(Float_t _eta, Float_t _pT) {

		if (!sfMap) return 1.;

		_pT = std::min(_pT, ptMax);
		_pT = std::max(_pT, ptMin);

		Int_t bin = sfMap->FindBin(_eta, _pT);
		return sfMap->GetBinContent(bin);
	};
};



struct binCenterTracker {
	//// http://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf


	/// key=bin up edge, value = (sum W*X, sum W*X^2, sum W, sum W^2)
	std::map<Double_t, std::tuple<Double_t, Double_t, Double_t, Double_t>> bins;
	Double_t 	minEdge;
	Double_t 	maxEdge;
	Bool_t 		excOverFlow = 1;
	Bool_t 		excUnderFlow = 1;

	binCenterTracker() {};

	binCenterTracker(std::vector<Double_t> _bins, Bool_t _excOverFlow, Bool_t _excUnderFlow) {
		init( _bins, _excOverFlow, _excUnderFlow);
	};

	void init(std::vector<Double_t> _bins, Bool_t _excOverFlow = 1, Bool_t _excUnderFlow = 1) {

		if (!std::is_sorted(_bins.begin(),_bins.end())) {
			std::cout<<"Error! bins are not sorted in ascending order."<<std::endl;
			return;
		}

		excOverFlow = _excOverFlow;
		excUnderFlow = _excUnderFlow;

		minEdge = _bins[0];
		maxEdge = _bins.back();

		for (UInt_t iEdge = 1; iEdge < _bins.size(); iEdge++) {
			bins[_bins[iEdge]] = std::make_tuple(0., 0., 0., 0.);
		}

	};

	void Fill(Double_t _xVal, Double_t _weight) {

		std::map<Double_t, std::tuple<Double_t, Double_t, Double_t, Double_t>>::iterator iBin;

		if (_xVal >= maxEdge) {

			if (excOverFlow) return;
			else iBin = std::prev(bins.end());

		}  else if (_xVal < minEdge) {

			if (excUnderFlow) return;
			else iBin = bins.begin();

		} else {

			iBin = bins.lower_bound(_xVal);

		}

		std::get<0>(iBin->second) += _xVal * _weight;
		std::get<1>(iBin->second) += _xVal*_xVal * _weight;
		std::get<2>(iBin->second) += _weight;
		std::get<3>(iBin->second) += _weight * _weight;
	};


	std::vector<std::pair<Double_t, Double_t>> getMeans() {

		std::vector<std::pair<Double_t, Double_t>> means;

		for (std::map<Double_t, std::tuple<Double_t, Double_t, Double_t, Double_t>>::iterator iBin = bins.begin(); iBin != bins.end(); iBin++) {

			if (std::get<2>(iBin->second) < std::numeric_limits<Double_t>::min()) {
				means.push_back(std::make_pair(std::numeric_limits<Double_t>::max() - 1000, std::numeric_limits<Double_t>::max() - 1000));
				continue;
			}

			Double_t mean = std::get<0>(iBin->second) / std::get<2>(iBin->second);

			Double_t unbiasedVariance = std::get<1>(iBin->second) / std::get<2>(iBin->second) - mean*mean; //// <X^2> - <X>^2
			unbiasedVariance = unbiasedVariance / (1. - std::get<3>(iBin->second)/std::pow(std::get<2>(iBin->second), 2)); /// bias correction : (<X^2> - <X>^2) / (1 - sum W^2 / (sum W)^2)

			Double_t nEff =	std::pow(std::get<2>(iBin->second),2)/std::get<3>(iBin->second);

			Double_t stdErr = std::sqrt(unbiasedVariance/nEff);

			means.push_back(std::make_pair(mean, stdErr));

			// std::cout<<"\t\t"<<std::get<0>(iBin->second)<<"\t"<<std::get<1>(iBin->second)<<"\t"<<std::get<2>(iBin->second)<<"\t"<<std::get<3>(iBin->second)<<std::endl;
		}

		return means;

	};

};



struct QuantileCorrector {
	QuantileCorrector () {};
	QuantileCorrector (std::string filePath, std::string graphName) {
		init(filePath, graphName);
	};

	void init(std::string filePath, std::string graphName) {
		delete correctionGraph;
		delete invertedGraph;
		correctionGraph = (TGraph*) getObjectFromFile(graphName, filePath);
		if (correctionGraph) std::cout<<"Loaded BDT corrections..."<<std::endl;
		else std::cout<<"Failed to load BDT corrections from "<<filePath<<" "<< graphName<<std::endl;
		createInvertedCorrections();
	};

	void createInvertedCorrections() {
		invertedGraph = new TGraph(correctionGraph->GetN(), correctionGraph->GetY(), correctionGraph->GetX());
	}

	Float_t getCorrectedBDTscore(Float_t uncorrectedBDTscore) {
		return correctionGraph->Eval(uncorrectedBDTscore);
	};

	Float_t getInverseCorrectedBDTscore(Float_t uncorrectedBDTscore) {
		return invertedGraph->Eval(uncorrectedBDTscore);
	};

	TGraph* correctionGraph = nullptr;
	TGraph* invertedGraph = nullptr;

	~QuantileCorrector () {
		delete correctionGraph;
		delete invertedGraph;
		correctionGraph = nullptr;
		invertedGraph = nullptr;
	};
};



class PileupReWeighting {
public:

	PileupReWeighting( ) { };

	PileupReWeighting( std::string mcFile, std::string dataFile, std::string mcHistName, std::string dataHistName, Bool_t _verbose=0) {
		init(mcFile, dataFile, mcHistName, dataHistName, _verbose);
	};

	Bool_t failedLoadingHists = 0;

	void init( std::string mcFile, std::string dataFile, std::string mcHistName, std::string dataHistName, Bool_t _verbose=0) {

		std::cout<<"\tCreating Pileup Reweighter with "<<std::endl<<
		"\t\t\tdata histogram "<< dataHistName<<" from file "<<dataFile<<std::endl<<
		"\t\t\tmc histogram "<< mcHistName<<" from file "<<mcFile<<std::endl;

		weights_ = (TH1*) getHistFromFile(dataHistName, dataFile);
		MC_distr_ = (TH1*) getHistFromFile(mcHistName, mcFile);

		if (!MC_distr_ || !weights_) {
			std::cout<<"Warning! PU reweighting histogram(s) not found! Will set PU weights to 1."<<std::endl;
			delete weights_;
			delete MC_distr_;
			failedLoadingHists = 1;
			return;
		}

		weights_->SetName(randStr().c_str());
		MC_distr_->SetName(randStr().c_str());

		int NBins = weights_->GetNbinsX();

		weights_->Scale( 1000000.0/ weights_->Integral(0, NBins+1));
		MC_distr_->Scale( 1000000.0/ MC_distr_->Integral(0, NBins+1));

		weights_->SetName(randStr().c_str());
		weights_->Divide(MC_distr_);
		delete MC_distr_;
		MC_distr_ = nullptr;

		if (_verbose) std::cout << "\t\tPileup weights: " << std::endl;

		for (int ibin = 1; ibin<NBins+1; ++ibin) {
			if (_verbose) std::cout << "\t\t\t" << ibin-1 << "\t" << weights_->GetBinContent(ibin) << std::endl;
		}
		std::cout<<"\tPileup weights calculated!"<<std::endl;
	};

	Double_t weight( Float_t n_int ) {
		if (failedLoadingHists) return 1.;
		Int_t bin = weights_->GetXaxis()->FindBin(n_int);
		return weights_->GetBinContent(bin);
	}

	~PileupReWeighting() {
		delete weights_;
		delete MC_distr_;
	};

protected:
	TH1* weights_	= nullptr;
	TH1* MC_distr_	= nullptr;
};


Double_t deltaR(Double_t _eta1, Double_t _phi1, Double_t _eta2, Double_t _phi2) {
	Double_t _deltaEta = std::abs(_eta1 - _eta2);
	Double_t _deltaPhi = deltaPhi(_phi1, _phi2);
	Double_t _deltaR = std::sqrt(_deltaEta*_deltaEta + _deltaPhi*_deltaPhi);
	return _deltaR;
};


Double_t deltaPhi(Double_t phi1, Double_t phi2) {
	Double_t tmpDeltaPhi = std::abs(phi2-phi1);
	Double_t minDeltaPhi = (tmpDeltaPhi > TMath::Pi()) ? (2*TMath::Pi() - tmpDeltaPhi) : tmpDeltaPhi;

	// Double_t minDeltaPhi = phi1 - phi2;
	// if (minDeltaPhi > TMath::Pi()) minDeltaPhi -= 2.*TMath::Pi();
	// else if (minDeltaPhi < -TMath::Pi()) minDeltaPhi += 2.*TMath::Pi();

	return minDeltaPhi;
};


std::string removeNonAlpha(std::string word) {
	word.erase(std::remove_if(word.begin(), word.end(),
		[](char ch) {
			return !::iswalnum(ch);
		}), word.end());
	return word;
};

std::string removeNonAlphaSmart(std::string word, std::string texts2eraseCSV) {

	word.erase(std::remove_if(word.begin(), word.end(),
		[](char ch) {
			Bool_t doRemove =(!std::iswalnum(ch)) && (ch != '.') && (ch != '_');
			return doRemove;
		}), word.end());

	word = findAndReplaceAll(word, ".", "p");

	std::vector<std::string> strings2erase 	=	split_string(texts2eraseCSV);

	for (std::string & iErase : strings2erase) {
		word = findAndReplaceAll(word, iErase, "");
	}

	return word;
};


template <class any_number>
std::string removeTrailingZeros(any_number number) {
	std::string str = std::to_string (number);
	str.erase( str.find_last_not_of('0') + 1, std::string::npos);
	str.erase(str.find_last_not_of('.') + 1, std::string::npos);
	if (str.length()>0 && !str.substr(str.length()-1,1).compare(".")) str.pop_back();
	return str;
};


Bool_t file_exists(std::string fileName) {
	// std::ifstream infile(fileName);
	// return infile.good();
	if (isDirectory(fileName)==1) return 0;
	struct stat buffer;
	return (stat (fileName.c_str(), &buffer) == 0);
};


Int_t mkdir(std::string dir_path) {
	if (trim_copy(dir_path).empty()) {
		std::cout<<"Empty path given! Cannot mkdir"<<std::endl;
		return -1;
	}
	std::string command = "mkdir -p " + dir_path;
	const int dir_err = system(command.c_str());

	// if(checkIfDirectory(dir_path)){
	// 	std::cout<<"Directory "<<dir_path<<" already exists"<<std::endl;
	// 	return 1;
	// }

	if (-1 == dir_err) {
		printf("Error creating directory!");
	} else std::cout <<"Created directory: " << dir_path << std::endl;
	return dir_err;
};


std::vector<std::string> getObjectList(std::string filepath, std::string objtype, std::vector<std::string> exclusion_list) {
	std::cout<<"Making object list..."<<std::endl;
	std::vector<std::string> obj_names;
	TFile rootfile(filepath.c_str(), "READ");
	TIter next(rootfile.GetListOfKeys());
	TKey *key;
	while ((key = (TKey*) next())) {
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom(objtype.c_str())) continue;
		TObject * g = key->ReadObj();
		std::string name = g->GetName();
		Bool_t to_exclude = (std::find(exclusion_list.begin(), exclusion_list.end(), name) != exclusion_list.end());
		if (to_exclude) {
			std::cout << "\t " << "Excluding "<< objtype<<" : "<< name << std::endl;
			continue;
		}
		std::cout << "\t " << "Added "<< objtype<<" : "<< name << std::endl;
		obj_names.push_back(name);
		g->Delete();
	}
	rootfile.Close();
	std::sort( obj_names.begin(), obj_names.end() );
	obj_names.erase( std::unique( obj_names.begin(), obj_names.end() ), obj_names.end());
	return obj_names;
};


Bool_t match(std::string _needle, std::string _haystack) {
	char const *needle = _needle.c_str();
	char const *haystack = _haystack.c_str();
	for (; *needle != '\0'; ++needle) {
		switch (*needle) {
		case '?': {
			if (*haystack == '\0')	return false;
			++haystack;
			break;
		}
	case '*': {
		if (needle[1] == '\0')	return true;
		size_t max = strlen(haystack);
		for (size_t i = 0; i < max; i++)
			if (match(needle + 1, haystack + i)) return true;
		return false;
	}
default:
	if (*haystack != *needle) return false;
	++haystack;
}
}
return *haystack == '\0';
};


Bool_t matchManyNeedles(std::vector<std::string> _needles, std::string _haystack) {
	Bool_t foundAneedle = 0;
	for (std::string & iNeedle : _needles) {
		if (match(iNeedle, _haystack)) {
			foundAneedle = 1;
			break;
		}
	}
	return foundAneedle;
};

std::string ReadNthLine(std::string filename, int N) {
	std::ifstream in(filename.c_str());
	std::string s;
	//for performance
	s.reserve(200);
	//skip N lines
	for (int i = 0; i < N; ++i) {
		std::getline(in, s);
	}
	std::getline(in,s);
	return s;
};


UInt_t countLines(std::string filename) {
	std::ifstream myfile(filename);
	// new lines will be skipped unless we stop it from happening:
	myfile.unsetf(std::ios_base::skipws);
	// count the newlines with an algorithm specialized for counting:
	UInt_t line_count = std::count(std::istream_iterator<char>(myfile),	std::istream_iterator<char>(), '\n');
	return line_count;
};


// std::vector<std::string> split_string(std::string _string, std::string _delimiter){
// 	std::vector<string> _split_string;
// 	boost::split(_split_string,_string,boost::is_any_of(_delimiter));
// 	return _split_string;
// };


std::vector<std::string> split_string(std::string _string, std::string _delimiter, Bool_t _trim) {
	size_t pos = 0;
	std::string token;
	std::vector<std::string> res;
	while ((pos = _string.find(_delimiter)) != std::string::npos) {
		token = _string.substr(0, pos);
		if (_trim) trim(token);
		_string.erase(0, pos + _delimiter.length());
		res.push_back(token);
	}
	if (_trim) trim(_string);
	res.push_back(_string);
	return res;
};


std::string get_cell(std::string filename, UInt_t row, UInt_t column, std::string _delimiter) {
	return split_string(ReadNthLine(filename, row), _delimiter)[column];
};


std::vector<std::string> splitpath(const std::string& str, const std::set<char> delimiters) {
	std::vector<std::string> result;

	char const* pch = str.c_str();
	char const* start = pch;
	for (; *pch; ++pch) {
		if (delimiters.find(*pch) != delimiters.end()) {
			if (start != pch) {
				std::string str(start, pch);
				result.push_back(str);
			} else {
				result.push_back("");
			}
			start = pch + 1;
		}
	}
	result.push_back(start);

	return result;
}


std::string getFileName(std::string _filepath) {
	constexpr char sep = '/';
	size_t i = _filepath.rfind(sep, _filepath.length());
	if (i != string::npos) {
		return (_filepath.substr(i+1, _filepath.length() - i));
	}
	return (_filepath);
};


std::string getDirPath(std::string _somePath) {
	std::string directory;
	const size_t last_slash_idx = _somePath.rfind('/');
	if (std::string::npos != last_slash_idx) {
		directory = _somePath.substr(0, last_slash_idx);
	}
	return directory;
};


std::vector<std::string> getNonemptyLines(std::string filepath) {
	std::ifstream infile(filepath);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(infile, line)) {
		trim(line);
		if (!line.empty()) lines.push_back(line);
	}
	return lines;
};


std::vector<std::string> getNonemptyLinesWithFilenameKeyword(std::string filepath, std::string keyword, std::string exclude) {
	std::ifstream infile(filepath);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(infile, line)) {
		std::string filename = getFileName(line);
		if (filename.empty()) continue;
		if (!match(keyword.c_str(), filename.c_str())) continue;
		if (!exclude.empty() && match(exclude.c_str(), filename.c_str())) continue;
		lines.push_back(line);
	}
	return lines;
};


std::vector<std::string> getLinesRegex(std::string _filepath, std::string _regexStr) {
	regex _regex(_regexStr);
	std::ifstream infile(_filepath);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(infile, line)) {
		if (!std::regex_match(line, _regex)) continue;
		lines.push_back(line);
	}
	return lines;
};


TH1* getHistFromFile(std::string _histName, std::string _filename, Bool_t _verbose, TDirectory *_dir) {

	if (!file_exists(_filename)) {
		std::cout<<"Error! File "<<_filename<<" does not exist!"<<std::endl;
		return nullptr;
	}

	TFile _file(_filename.c_str(), "READ");
	if (!objectExists(&_file, _histName)) {
		_file.Close();
		std::cout<<"Error! Histogram "<<_histName<<" not found in file "<<_filename<<std::endl;
		return nullptr;
	}

	TH1* _hist = (TH1*) _file.Get(_histName.c_str());
	if (!_dir)_hist->SetDirectory(0);
	else _hist->SetDirectory(_dir);
	_file.Close();
	if (_verbose)std::cout<<"Loaded histogram "<<_histName<<" (N="<<_hist->GetEntries() <<") from file "<<_filename<<std::endl;
	if (!_hist) std::cout<<"Error! Histogram "<< _histName<<" not found in file "<<_filename<<std::endl;
	return _hist;
};


TObject *getObjectFromFile(std::string _objectName, std::string _filename) {

	if (!file_exists(_filename)) {
		std::cout<<"Error! File "<<_filename<<" not found!"<<std::endl;
		return nullptr;
	}

	TFile _file(_filename.c_str(), "READ");
	TObject* _object = (TObject*) _file.Get(_objectName.c_str());
	// _object->SetDirectory(0);
	_file.Close();
	if (!_object) {
		std::cout<<"Error reading "<<_objectName<< " from file "<<_filename<<std::endl;
		return nullptr;
	}
	std::cout<<"\t\tLoaded TObject "<<_objectName<<" from file "<<_filename<<std::endl;
	return _object;
};


TH1* rebinHist(TH1* _hist, Double_t _statUnc) {
	std::vector<Double_t> _goodBins = getGoodBins(_hist, _statUnc);
	if (_goodBins.size()<2) return _hist;
	std::string _newname = "rebinned_" + (std::string)_hist->GetName();
	TH1* _rebinnedHist = (TH1*) _hist->Rebin(_goodBins.size()-1, _newname.c_str(), _goodBins.data());
	std::cout<<"\t\t\tRebinned "<<_hist->GetName()<<". New bins:";
	for (Double_t iBin : _goodBins) {
		std::cout<<"\t"<<iBin;
	}
	std::cout<<std::endl;
	// _rebinnedHist->Scale(1.,"width");
	return _rebinnedHist;
};


TH1* rebinHist(TH1* _hist, std::vector<Double_t> _newBins) {
	std::string _newname = "rebinned_" + (std::string)_hist->GetName();
	TH1 *_hist_rebinned = (TH1*) _hist->Rebin(_newBins.size()-1, _newname.c_str(), _newBins.data());
	std::cout<<"\t\t\tRebinned "<<_hist->GetName()<<". New bins:";

	for (Double_t iBin : _newBins) {
		std::cout<<"\t"<<iBin;
	}

	_hist->Delete();
	return _hist_rebinned;
};


TH1* rebinHist(TH1* _hist, TH1* _templateHist) {
	TH1 *_hist_rebinned = (TH1*) _templateHist->Clone((std::string(_hist->GetName())+"_rebinned").c_str());
	_hist_rebinned->Reset();
	for (Int_t i = 1; i <= _hist->GetNbinsX(); i++) {
		Int_t iNewBin = _hist_rebinned->GetXaxis()->FindBin(_hist->GetXaxis()->GetBinCenter(i));
		Double_t rebinContent = _hist_rebinned->GetBinContent(iNewBin) + _hist->GetBinContent(i); 		//// Add existing content in new
		_hist_rebinned->SetBinContent(iNewBin, rebinContent);
	}

	return _hist_rebinned;
};


void removeBinErrors(TH1* _hist){
	for(Int_t iX = 1; iX <= _hist->GetNbinsX(); iX++){
		for(Int_t iY = 1; iY <= _hist->GetNbinsX(); iY++){
			_hist->SetBinError(iX, iY, 0.);
			
		}	
	}
};


TH1* rebinNHist(TH1* _hist, Int_t _N) {
	std::string _newname = "rebinned_" + std::to_string(_N) + "_" + (std::string)_hist->GetName();
	TH1 *_hist_rebinned = (TH1*) _hist->Rebin(_N, _newname.c_str());
	_hist->Delete();
	return _hist_rebinned;
};

Double_t sumNextNbins(TH1* _hist, Int_t _n, Int_t _curr) {
	Double_t _NbinSum = 0.;
	Int_t _countEnd = (_hist->GetNbinsX() <= (_curr + _n)) ? _hist->GetNbinsX() : _curr + _n;
	for (Int_t i = _curr; i < _countEnd+1; i++) {
		_NbinSum += _hist->GetBinContent(i);
	}
	return _NbinSum;
};


std::vector<Double_t> getGoodBins(TH1* _hist, Double_t _statUnc, Double_t _reScale, Int_t _nbinPar) {

	TH1 *reScaled = (TH1*)_hist->Clone("reScaled");
	if (_reScale > 0.) {
		reScaled->Scale(_reScale / reScaled->Integral(0, reScaled->GetNbinsX()));
		_hist = reScaled;
	}
	UInt_t _nBins = _hist->GetXaxis()->GetNbins();
	if (_nBins == 0) {
		std::cout<<"\t Hist "<<_hist->GetName()<<" has no bins!"<<std::endl;
		return {};
	}
	_hist->Sumw2();
	std::vector<Double_t> good_bins;
	good_bins.push_back(_hist->GetXaxis()->GetBinLowEdge(1));
	Double_t _runningBinSum = 0.;
	UInt_t _first50pcBins = std::ceil(0.5 *_nBins);
	for (UInt_t i = 1; i < _nBins+1; i++) {
		Double_t _bincontent = _hist->GetBinContent(i);
		Double_t _nextNbincontent = sumNextNbins(_hist,_nbinPar,i);
		// if(_bincontent != _bincontent) return {};
		Double_t _binUpedge = _hist->GetXaxis()->GetBinUpEdge(i);
		Double_t _binLowedge = _hist->GetXaxis()->GetBinLowEdge(i);
		//(i <=_first50pcBins && _bincontent == 0.)
		// if((_nextNbincontent == 0. && i < _nBins) || (_bincontent == 0.)){
		// if((_nextNbincontent == 0. && i < _nBins)){
		// 	if(good_bins.size()>0){
		// 		Double_t _lastElement = good_bins.back();
		// 		if(_binLowedge -_lastElement > 0.000000001) good_bins.push_back(_binLowedge);
		// 	}
		// 	good_bins.push_back(_binUpedge);
		// 	_runningBinSum = 0.;
		// 	continue;
		// }
		_runningBinSum += _bincontent;
		if (1/std::sqrt(_runningBinSum) < _statUnc) {
			good_bins.push_back(_binUpedge);
			_runningBinSum =0;
		}
		if (i < _nBins) continue;
		Double_t _lastElement = good_bins.back();
		if (_binUpedge -_lastElement < 0.00000001) continue;
		good_bins.push_back(_binUpedge);
		_runningBinSum =0;
		if (1/std::sqrt(_runningBinSum) > _statUnc) continue;
		if (good_bins.size()==1) continue;
		//erase second last element
		good_bins.erase(good_bins.begin() + good_bins.size()-2);
	};
	reScaled->Delete();
	return good_bins;
};



Double_t getSumW(std::string _cutflowfile) {
	std::string sumWline = ReadNthLine(_cutflowfile, 3);
	std::string sumWstring = (split_string(sumWline, ",")).at(0);
	Double_t _sumW = std::stod(sumWstring);
	std::cout<<"\tFrom "<<_cutflowfile<<" sumW="<<_sumW<<std::endl;
	return _sumW;
};


// trim from start (in place)
void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		return !std::isspace(ch);
	}));
}


// trim from end (in place)
void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		return !std::isspace(ch);
	}).base(), s.end());
}


// trim from both ends (in place)
void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}


// trim from start (copying)
std::string ltrim_copy(std::string s) {
	ltrim(s);
	return s;
}


// trim from end (copying)
std::string rtrim_copy(std::string s) {
	rtrim(s);
	return s;
}


// trim from both ends (copying)
std::string trim_copy(std::string s) {
	trim(s);
	return s;
};


Double_t ams(Double_t _s, Double_t _b) {
	return std::sqrt(2*(1+_s+_b)*std::log(1+_s/_b)-2*_s);
};


template <typename T>
std::string to_string_with_precision(const T a_value, const int n) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(n) << a_value;
	return out.str();
};


std::string getUnit(TH1* _hist) {
	std::string ystring = _hist->GetYaxis()->GetTitle();
	return getUnit(ystring);
};


std::string getUnit(std::string ystring) {
	size_t lastSlash = ystring.find_last_of("/");
	if (lastSlash == std::string::npos) return "";
	std::string _unit = ystring.substr(lastSlash+1);
	trim(_unit);
	std::string number = first_numberstring(_unit);
	std::string unit_text = _unit;
	boost::replace_all(unit_text, number, "");
	trim(unit_text);
	if (unit_text.empty()) return std::to_string(1);
	else return unit_text;
};


std::string eraseUnit(std::string ystring) {
	size_t lastSlash = ystring.find_last_of("/");
	if (lastSlash == std::string::npos) return ystring;
	std::string _unit = ystring.substr(lastSlash);
	boost::replace_all(ystring, _unit, "");
	return ystring;
}


std::string first_numberstring(std::string const & str) {
	std::size_t const n = str.find_first_of("0123456789.");
	if (n != std::string::npos) {
		std::size_t const m = str.find_first_not_of("0123456789.", n);
		return str.substr(n, m != std::string::npos ? m-n : m);
	}
	return std::string();
};


void closeTChain(TChain * & _chain) {
	if (!_chain) {
		std::cout<<"Error! TChain is null!"<<std::endl;
	}
	TFile *file = _chain->GetCurrentFile();
	_chain->SetDirectory(0);
	if (file) delete file;
	_chain = nullptr;
};

void closeTTree(TTree * & _chain) {
	if (!_chain) {
		std::cout<<"Error! TChain is null!"<<std::endl;
	}
	TFile *file = _chain->GetCurrentFile();
	_chain->SetDirectory(0);
	if (file) file->Close();
	_chain = nullptr;
};


std::vector<Double_t> getXbins(TH1* hist) {
	UInt_t nbins = hist->GetNbinsX();
	if (nbins == 0) {
		std::cout << "Error : No Bins!" << std::endl;
		return {};
	};
	std::vector<Double_t> bins;
	TAxis* axis = hist->GetXaxis();
	bins.push_back(axis->GetBinLowEdge(1));
	for (UInt_t i = 1; i < nbins + 1; i++) {
		bins.push_back(axis->GetBinUpEdge(i));
	}
	return bins;
};


void setFrameColor(TAxis* _axis, std::string _color) {
	Int_t color_val = TColor::GetColor(_color.c_str());
	_axis->SetAxisColor(color_val);
	_axis->SetLabelColor(color_val);
	_axis->SetTitleColor(color_val);
};


void setFrameColor(TH1* _hist, std::string _color) {
	setFrameColor(_hist->GetXaxis(), _color);
	setFrameColor(_hist->GetYaxis(), _color);
};


void setFrameColor(THStack* _stack, std::string _color) {
	setFrameColor(_stack->GetXaxis(), _color);
	setFrameColor(_stack->GetYaxis(), _color);
};


TTree *loadTreeFromFile(std::string _treeName, std::string _filename) {
	TFile *_file = new TFile(_filename.c_str(), "READ");
	if (!_file) {
		std::cout<<"Error! Cannot read file "<<_filename<<std::endl;
		return nullptr;
	}

	TTree *_tree = (TTree*) _file->Get(_treeName.c_str());
	if (!_tree) {
		std::cout<<"Error! Cannot read TTree "<<_file<<" from successfully loaded TFile "<<_filename <<std::endl;
		return nullptr;
	}

	std::cout<<"Successfully loaded TTree "<<_treeName<<" from TFile "<<_filename <<std::endl;
	std::cout<<"Remember to call TTree*->GetDirectory()->Close() once the tree is not needed any more!" <<std::endl;

	return _tree;
};


Char_t isDirectory(std::string filePath, Bool_t _verbose) {
	DIR* dir = opendir(filePath.c_str());
	if (dir) {
		closedir(dir);
		return 1;
	} else if (ENOENT == errno) {
		return 0;
	} else {
		if (_verbose) std::cout<<"Error checking path : "<<filePath<<std::endl;
		return -1;
	}
}


TGraph* removeErrors(TGraphErrors* graphToConvert) {

	TGraph *noErrorGraph = new TGraph(graphToConvert->GetN(), graphToConvert->GetX(), graphToConvert->GetY());

	return noErrorGraph;

};


void addPointToGraph(TGraphAsymmErrors * _graph, Double_t _x, Double_t _xErrL, Double_t _xErrH, Double_t _y, Double_t _yErrL, Double_t _yErrH) {
	UInt_t nExistingPts = _graph->GetN();
	_graph->SetPoint(nExistingPts, _x, _y);
	_graph->SetPointError(nExistingPts, _xErrL, _xErrH, _yErrL, _yErrH);
}

void addPointToGraph(TGraph * _graph, Double_t _x, Double_t _y) {
	_graph->SetPoint(_graph->GetN(), _x, _y);
};


void addPointToGraph(TGraphErrors * _graph, Double_t _x, Double_t _ex, Double_t _y, Double_t _ey) {
	Int_t lastGraphPoint = _graph->GetN();
	_graph->SetPoint(lastGraphPoint, _x, _y);
	_graph->SetPointError(lastGraphPoint, _ex, _ey);
};

void graphStats(TGraphAsymmErrors* graph, Double_t &mean, Double_t &stdev) {
	UInt_t Npoints = graph->GetN();
	std::vector<Double_t> xVals(graph->GetY(), graph->GetY() + Npoints), xErrsL(graph->GetEYlow(), graph->GetEYlow() + Npoints), xErrsH(graph->GetEYhigh(), graph->GetEYhigh() + Npoints);
	Double_t W, sumW = 0, sumW2 = 0, sumXW = 0, sumDev2W = 0;
	for (UInt_t i = 0; i < Npoints; i++) {
		if (xErrsL[i] * xErrsH[i] > 0) {
			W = 1 / (xErrsL[i] * xErrsH[i]);
			sumW += W;
			sumXW += xVals[i] * W;
		}
	}
	if (sumW > 0) {
		mean = sumXW / sumW;
		for (UInt_t i = 0; i < Npoints; i++) {
			if (xErrsL[i] * xErrsH[i] > 0) {
				W = 1 / (xErrsL[i] * xErrsH[i]);
				sumW2 += W * W;
				sumDev2W += (xVals[i] - mean) * (xVals[i] - mean) * W;
			}
		}
		stdev = TMath::Sqrt(sumDev2W / (sumW - sumW2 / sumW));
	} else {
		mean = -999;
		stdev = -999;
	}
};


TH1D* graph2hist(TGraphAsymmErrors* graph, UInt_t ndivs, Double_t ylow, Double_t yhigh) {
	std::string name = graph->GetName();
	name.append("_hist_" + to_string_with_precision(rand() % 1000, 0));
	UInt_t Npoints = graph->GetN();

	std::vector<Double_t> xVals(graph->GetY(), graph->GetY() + Npoints), xErrsL(graph->GetEYlow(), graph->GetEYlow() + Npoints), xErrsH(graph->GetEYhigh(), graph->GetEYhigh() + Npoints);

	TH1D* hist = new TH1D("", "", ndivs, ylow, yhigh);
	hist->SetFillColor(graph->GetLineColor());
	hist->SetLineColor(graph->GetLineColor());
	hist->SetMarkerColor(graph->GetLineColor());

	Double_t w = 1;
	for (UInt_t i = 0; i < Npoints; i++) {
		if (!(xErrsL[i] <= 0  || xErrsH[i] <= 0))w = 1 / (xErrsL[i] * xErrsH[i]);
		else continue;
		hist->Fill(xVals[i], w);
	}
	Npoints = hist->GetEffectiveEntries();
	hist->Scale(Npoints / hist->Integral());
	return hist;
};


Int_t writeToFile(TObject *_object, std::string _filePath, std::string _mode, Bool_t _verbose, std::string _writeName, std::string _tdir) {

	std::string writeDir = getDirPath(_filePath);
	if (!isDirectory(writeDir)) mkdir(writeDir);

	TFile _outfile(_filePath.c_str(), _mode.c_str());

	if (_outfile.IsZombie()) {
		std::cout<<"\t\tFailed to open TFile "<<_filePath<<" in the mode "<<_mode<<std::endl;
		return -1;
	}

	TDirectory* theDirectory = nullptr;

	if (_object) {
		_outfile.cd();
		if (!_tdir.empty()) {
			theDirectory = _outfile.GetDirectory(_tdir.c_str());
			if (theDirectory == nullptr) theDirectory = _outfile.mkdir(_tdir.c_str());
			theDirectory->cd();
		}
		if (_writeName.empty()) {
			_object->Write(_object->GetName());
			if (_verbose) std::cout<<"\t\tWritten "<< _object->GetName()<<" to file "<<_filePath<<std::endl;
		} else {
			_object->Write(_writeName.c_str());
			if (_verbose) std::cout<<"\t\tWritten "<< _writeName<<" to file "<<_filePath<<std::endl;
		}

	}
	_outfile.Close();
	delete theDirectory;
	return 0;
};

Int_t writeHistToFile(TH1 *_object, std::string _filePath, std::string _mode, Bool_t _verbose, std::string _writeName) {

	TH1* objCpy = (TH1*) _object->Clone(_object->GetName());
	objCpy->ResetAttFill();
	objCpy->ResetAttLine();
	objCpy->ResetAttMarker();
	Int_t retVal = writeToFile(objCpy, _filePath, _mode, _verbose, _writeName);
	delete objCpy;
	return retVal;
};


std::string findAndReplaceAll(std::string data, std::string toSearch, std::string replaceStr) {
	// Get the first occurrence
	size_t pos = data.find(toSearch);

	// Repeat till end is reached
	while ( pos != std::string::npos) {
		// Replace this occurrence of Sub String
		data.replace(pos, toSearch.size(), replaceStr);
		// Get the next occurrence from the current position
		pos =data.find(toSearch, pos + replaceStr.size());
	}

	return data;
};

TChain *openTree(std::string _inFile, std::string _treeName, Bool_t _verbose) {
	TChain* theTree = openTChain(std::vector<std::string> {_inFile}, _treeName, _verbose);
	return theTree;
};

TChain *openTChain(std::string _chainListFile, std::string _treeName, Bool_t _verbose) {
	if (!file_exists(_chainListFile)) {
		std::cout<<"Error! Chain list file does not exist:"<<_chainListFile<<std::endl;
		exit(EXIT_FAILURE);
	}
	if (_verbose) std::cout<<"\tMaking TChain with root files listed in "<<_chainListFile<<std::endl;
	std::vector<std::string> _ntuples = getNonemptyLines(_chainListFile);
	if (_ntuples.size() == 0) {
		std::cout<<"Error! File list "<<_chainListFile<<" is invalid!"<<std::endl;
		return nullptr;
	}

	if (_ntuples.empty()) {
		std::cout<<"Error! TChain list is empty!"<<std::endl;
		return nullptr;
	}

	return openTChain(_ntuples, _treeName, _verbose);
};


TChain *openTChain(std::vector<std::string> _chainList, std::string _treeName, Bool_t _verbose) {

	if (_chainList.empty()) {
		std::cout<<"Error! TChain list is empty!"<<std::endl;
		return nullptr;
	}

	if (_treeName.empty()) {
		TFile _testF(_chainList[0].c_str(), "READ");
		if (!_testF.GetListOfKeys()) {
			std::cout<<"Error! Cannot open TChain. No key found in file "<<_chainList[0]<<std::endl;
			exit(EXIT_FAILURE);
		}
		TIter lNext(_testF.GetListOfKeys()) ;
		TObject* lObj ;
		Char_t treeCounter = 0;
		while ((lObj = (TObject*)lNext())) {
			if (treeCounter > 1) {
				std::cout<<"Error! No TTree name provided & there are >1 trees in the file "<<_chainList[0]<<std::endl;
				return nullptr;
			}
			if (lObj->InheritsFrom(TTree::Class())) {
				_treeName = lObj->GetName();
				treeCounter++;
			}
		}
		_testF.Close();
	}

	if (_treeName.empty()) {
		std::cout<<"Error! No TTree name found!"<<std::endl;
		return nullptr;
	}

	TChain *_bChain = new TChain(_treeName.c_str());
	if (_verbose) std::cout<<"\tMaking TChain from trees "<<_treeName<<std::endl;
	for (auto & _ntuple : _chainList) {
		if (!file_exists(_ntuple)) {
			std::cout<<"Error! File does not exist "<<_ntuple<<std::endl;
			closeTChain(_bChain);
			return nullptr;
		};
		_bChain->Add(_ntuple.c_str());
		if (_verbose) std::cout << "\tAdded file "<< _ntuple <<std::endl;
	}
	if (_verbose) std::cout<<"\tBuilt TChain!"<<std::endl;

	return _bChain;
};


// TChain *openTChainWithFilesInDir(std::string _dirPath, std::string _treeName) {
// 	if (isDirectory(_dirPath)!=1) {
// 		std::cout<<"Error! Directory "<<_dirPath<<" doesn't exist"<<std::endl;
// 		return nullptr;
// 	}

// 	std::cout<<"Making TChain from trees "<<_treeName<<" in directory "<<_dirPath<<std::endl;
// 	std::vector<std::string> _fileList = listFilesInDir(_dirPath, ".*\.root", 1);
// 	// TChain *_bChain = openTChain(_fileList, _treeName);
// 	return openTChain(_fileList, _treeName);
// };


std::vector<std::pair<std::string, std::string>> getBranchList(std::string _treePath, std::string _treeName, Bool_t _verbose) {
	TChain *bChain = openTChain(std::vector<std::string> {_treePath}, _treeName);
	if (!bChain) {
		std::cout<<"Error! Could no create TChain!"<<std::endl;
		return {};
	}

	TObjArray *bArray = (TObjArray*)bChain->GetListOfBranches()->Clone();
	bArray->SetOwner(kFALSE);
	// bArray->Sort();

	std::vector<std::pair<std::string, std::string>> bList;

	std::cout<<"List of branches in TTree "<<_treeName<<" in file "<<_treePath<<":"<<std::endl;
	TIter iBr(bArray);
	TObject* bObj = nullptr;
	while ((bObj = (TObject*)iBr())) {
		std::string bName = bObj->GetName();
		std::string bType = bChain->GetLeaf(bName.c_str())->GetTypeName();
		bList.emplace_back(bType, bName);
		if (_verbose) std::cout<<"\t"<<bType<<",\t"<<bName<<std::endl;
	}

	closeTChain(bChain);

	return bList;
};


std::vector<std::string> listFilesInDir(std::string _dirPath, std::string _regexStr, Bool_t _verb) {

	if (!isDirectory(_dirPath)) {
		std::cout<<"Error! Directory ("<<_dirPath<<") does not exist!"<<std::endl;
		return {};
	}

	std::vector<std::string> _fileList;

	DIR           *dirp;
	struct dirent *directory;

	regex regex_;

	if (!_regexStr.empty()) {
		regex regextmp_(_regexStr);
		regex_ = regextmp_;
	}

	if (_verb) std::cout <<"Files in directory "<< _dirPath<<" matching  \""<<_regexStr<< "\" :" << std::endl;
	dirp = opendir(_dirPath.c_str());
	if (dirp) {
		while ((directory = readdir(dirp)) != NULL) {
			if (!_regexStr.empty()) {
				if (!std::regex_match(directory->d_name, regex_)) continue;
			}
			if (_verb) std::cout <<"\t\t"<< getFileName(directory->d_name) << std::endl;
			std::string _filePath = _dirPath + "/" + directory->d_name;
			_fileList.push_back(_filePath);
		}
		closedir(dirp);
	}

	return _fileList;
};


Bool_t matchRegex(std::string _str,std::string _regexStr) {
	regex regex_(_regexStr);
	return std::regex_match(_str, regex_);
};


std::string getTreeNameInFile(std::string _filePath) {
	std::string _treeName="";
	TFile _testF(_filePath.c_str(), "READ");
	if (!_testF.GetListOfKeys()) {
		std::cout<<"Error! Cannot open TChain. No key found in file "<<_filePath<<std::endl;
	}
	TIter lNext(_testF.GetListOfKeys()) ;
	TObject* lObj ;
	Char_t treeCounter = 0;
	while ((lObj = (TObject*)lNext())) {
		if (treeCounter > 1) {
			std::cout<<"Error! There are >1 trees in the file "<<_filePath<<std::endl;
			return "";
		}
		if (lObj->InheritsFrom(TTree::Class())) {
			_treeName = lObj->GetName();
			treeCounter++;
		}
	}
	_testF.Close();
	return _treeName;
};


// 	mergeBins
// 	Adds histograms for MC samples of the same process produced in different phase space bins
// 	_inFiles : list of root files containing histograms
// 	_histName : the histogram to be merged
// 	_xsecMap : CSV file listing samples and cross sections
// 	root file name is used to look up cross section
TH1F *mergeBins(std::vector<std::string> _inFiles, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol, Int_t _xSecCol, std::string _path) {

	TH1F *_mergedHist = (TH1F*) getHistFromFile(_histName, _path + _inFiles[0]);
	_mergedHist->Reset("ICESM");
	_mergedHist->SetDirectory(0);
	_mergedHist->Sumw2();
	std::string _newName = _mergedHist->GetName();
	_newName += "_mergedBins";
	_mergedHist->SetName(_newName.c_str());

	for (const auto & _file : _inFiles) {
		std::string _sampleName = getFileName(_file);
		_sampleName = findAndReplaceAll(_sampleName, ".root", "");

		TH1F *_sumWHist = (TH1F*) getHistFromFile(_sumWeightsHistname, _path + _file);
		Double_t _sumW = _sumWHist->GetBinContent(1);
		_sumWHist->Delete();

		Double_t _xSection = std::stod(vLookup(_sampleName, _xsecMap, _nameCol, _xSecCol));

		TH1F *_binHist = (TH1F*) getHistFromFile(_histName, _path + _file);
		_binHist->Scale(_xSection/_sumW);

		Double_t _integral = _binHist->Integral();

		_mergedHist->Add(_binHist);
		_binHist->Delete();

		std::cout<<"\t\t"<<_file<<":\t xSection = "<<_xSection<<"\tSumW = "<<_sumW<<"\tIntegral = "<< _integral<<std::endl;

	}

	std::cout<<"\t\t\tMerged all bins! Integral = "<<_mergedHist->Integral()<<std::endl;

	return _mergedHist;
};


// 	mergeBins
// 	Adds histograms for MC samples of the same process produced in different phase space bins
// 	_fileList : list of root files containing histograms
// 	_histName : the histogram to be merged
// 	_xsecMap : CSV file listing samples and cross sections
// 	root file name is used to look up cross section
TH1F *mergeBins(std::string _fileList, std::string _histName, std::string _sumWeightsHistname, std::string _xsecMap, Int_t _nameCol, Int_t _xSecCol, std::string _path) {
	std::cout<<"Merging TH1F "<<_histName<<" using file list "<< _fileList<<std::endl;
	std::vector<std::string> _inFiles;
	if (!isROOTfile(_fileList))_inFiles = getNonemptyLines(_fileList);
	else _inFiles = {_fileList};

	std::sort(_inFiles.begin(), _inFiles.end());

	return mergeBins(_inFiles, _histName, _sumWeightsHistname, _xsecMap, _nameCol, _xSecCol, _path);
};


std::string vLookup(std::string _lookupKey, std::string _inFile, Int_t _lookupCol, Int_t _valCol, Bool_t _regex, std::string _delimiter, Bool_t _silent) {
	CSVReader reader(_inFile, _delimiter);
	std::vector<std::vector<std::string>> data_matrix = reader.getData();

	Int_t _matchedRow = -999;
	for (UInt_t i = 0; i < data_matrix.size(); i++) {
		std::string _searchCell = data_matrix[i][_lookupCol];

		if (_regex) {
			if (matchRegex(_searchCell, _lookupKey)) {
				_matchedRow = i;
			}
		} else {
			if (_searchCell == _lookupKey) {
				_matchedRow = i;
			}
		}
		if (_matchedRow>-1) break;
	}
	if (_matchedRow < 0) {
		if (!_silent) std::cout<<"vLookup failed to "<<(_regex ? "match " : "find ")<<"key "<<_lookupKey<<" in file "<<_inFile<<std::endl;
		return "";
	}

	return data_matrix[_matchedRow][_valCol];
};


std::string hvLookup(std::string _rowLookupKey, std::string _colLookupKey, std::string _inFile, Int_t _lookupRow, Int_t _lookupCol, Bool_t _regex, std::string _delimiter, Bool_t _silent) {

	CSVReader reader(_inFile, _delimiter);

	std::vector<std::string> _lookupRowList = reader[_lookupRow];
	Int_t _matchedColumn = -999;
	for (UInt_t i = 0; i < _lookupRowList.size(); i++) {

		if (_regex) {
			if (matchRegex(_lookupRowList[i], _colLookupKey)) {
				_matchedColumn = i;
			}
		} else {
			if (_lookupRowList[i] == _colLookupKey) {
				_matchedColumn = i;
			}
		}
		if (_matchedColumn>-1) break;
	}

	_lookupRowList.clear();

	if (_matchedColumn < 0) {
		if (!_silent) std::cout<<"Error! Column lookup key \""<< _colLookupKey <<"\" not found in  file " <<_inFile<<std::endl;
		return "";
	}


	Int_t _matchedRow = -999;
	for (UInt_t i = 0; i < reader.size(); i++) {

		if (_regex) {
			if (matchRegex(reader[i][_lookupCol], _rowLookupKey)) {
				_matchedRow = i;
			}
		} else {
			if (reader[i][_lookupCol] == _rowLookupKey) {
				_matchedRow = i;
			}
		}
		if (_matchedRow>-1) break;
	}

	if (_matchedRow < 0) {
		if (!_silent) std::cout<<"Error! Row lookup key \""<< _rowLookupKey <<"\" not found in  file " <<_inFile<<std::endl;
		return "";
	}

	return reader[_matchedRow][_matchedColumn];
};


Bool_t isROOTfile(std::string _filepath) {
	if (!file_exists(_filepath)) {
		return 0;
	}
	if (_filepath.substr(_filepath.find_last_of(".") + 1) == "root") {
		if (TFile(_filepath.c_str()).GetListOfKeys() == nullptr) {
			return 0;
		}
		return 1;
	}
	return 0;
};


std::vector<Float_t> getXlimits(std::vector<TH1*> _hists, Float_t _binThreshold) {
	std::vector<Float_t> _limits;
	for (auto _hist: _hists) {
		Double_t first = _hist->GetXaxis()->GetBinLowEdge(_hist->FindFirstBinAbove(_binThreshold));
		Double_t last = _hist->GetXaxis()->GetBinUpEdge(_hist->FindLastBinAbove(_binThreshold));
		_limits.push_back(first);
		_limits.push_back(last);
	}

	std::vector<Float_t> _max_min;
	_max_min.push_back(*std::max_element(_limits.begin(), _limits.end()));
	_max_min.push_back(*std::min_element(_limits.begin(), _limits.end()));
	return _max_min;
};


void clearHeap() {
	TList* obList = gDirectory->GetList();
	for (auto obj : *obList) {
		obj->Delete();
	}
};


Double_t weightedYmean(TH1 *_hist) {
	Double_t _weightedSumY = 0.;
	Double_t _weightSum = 0.;
	for (Int_t i = 1; i <= _hist->GetNbinsX(); i++) {
		if (!(_hist->GetBinContent(i) > 0.)) continue;
		// if(!(_hist->GetBinError(i) > 0.)) continue;
		// Double_t _binWeight = 1./(_hist->GetBinError(i) * _hist->GetBinError(i));
		Double_t _binWeight = 1.;
		_weightSum += _binWeight;
		_weightedSumY += _binWeight * _hist->GetBinContent(i);
	}
	return (_weightedSumY/_weightSum);
};


Double_t weightedYspread(TH1 *_hist) {
	Double_t _weightedSumDeltaY2 = 0.;
	Double_t _weightSum = 0.;
	Double_t _weightedYmean = weightedYmean(_hist);

	for (Int_t i = 1; i <= _hist->GetNbinsX(); i++) {
		if (!(_hist->GetBinContent(i) > 0.)) continue;
		// if(!(_hist->GetBinError(i) > 0.)) continue;
		// Double_t _binWeight = 1./(_hist->GetBinError(i) * _hist->GetBinError(i));
		Double_t _binWeight = 1.;
		_weightSum += _binWeight;
		_weightedSumDeltaY2 += _binWeight * (_hist->GetBinContent(i) - _weightedYmean) * (_hist->GetBinContent(i) - _weightedYmean);
	}
	return (_weightedSumDeltaY2/_weightSum);
};


Bool_t branchExists(std::string _branchName, TTree *_tree) {
	if (!_tree) return 0;
	TBranch* br = (TBranch*) _tree->GetListOfBranches()->FindObject(_branchName.c_str());
	if (br)	return 1;
	else return 0;
};


Float_t getMean(std::vector<Float_t> _set) {
	if (_set.empty()) return -9999.;
	Double_t _sum = 0.;
	for (Float_t _num : _set) {
		_sum += _num;
	}
	return _sum/((Double_t)_set.size());
};


template <class ObjType>
ObjType copyObjectDeleteSrc(ObjType *_original) {
	ObjType _copy(*_original);
	_original->Delete();
	return _copy;
};


template<typename T1, typename T2>
Int_t findSecondaryIndex(T1 searchIndex, std::vector<T2> container) {
	T2 secondaryIndex = -999;
	for (UInt_t i = 0; i < container.size(); i++) {
		if ( container[i] ==  (T2) searchIndex) {
			secondaryIndex = i;
			break;
		}
	}
	return secondaryIndex;
};


Short_t findSecondaryIndex(Short_t searchIndex, std::vector<Short_t> container) {
	Short_t secondaryIndex = -999;
	for (UInt_t i = 0; i < container.size(); i++) {
		if ( container[i] ==  searchIndex) {
			secondaryIndex = i;
			break;
		}
	}
	return secondaryIndex;
};


std::string getCurrentTime() {
	std::chrono::time_point<std::chrono::system_clock> _now = std::chrono::system_clock::now();
	std::time_t _now_ = std::chrono::system_clock::to_time_t(_now);
	std::string c_time = std::ctime(&_now_);
	c_time.erase(std::remove(c_time.begin(), c_time.end(), '\n'), c_time.end());
	return c_time;
};


template <typename anytype>
void eraseElement(std::vector<anytype> & _vec, UInt_t _Index2Erase) {
	_vec.erase(_vec.begin() + _Index2Erase);
};


Double_t getCategoryBoundary(TH1 *_signal, TH1*_background) {

	std::cout<<"Signal integral:"<<_signal->Integral()<<std::endl<<"Background integral:"<<_background->Integral()<<std::endl;

	Int_t optimalBoundaryBin = -999;
	Int_t optimalJointSignificance = -999.;

	for (Int_t i = 1; i <= _signal ->GetNbinsX(); i++) {

		Double_t sIntegralBelow = _signal->Integral(0, i);
		Double_t bIntegralBelow = _background->Integral(0, i);

		Double_t sIntegralAbove = _signal->Integral(i+1, _signal ->GetNbinsX()+1);
		Double_t bIntegralAbove = _background->Integral(i+1, _signal ->GetNbinsX()+1);

		if (sIntegralAbove < 0.1*sIntegralBelow) continue;

		Double_t iJointSignificance = std::sqrt(std::pow(sIntegralBelow/std::sqrt(sIntegralBelow + bIntegralBelow), 2) + std::pow(sIntegralAbove/std::sqrt(sIntegralAbove + bIntegralAbove), 2));

		if (optimalJointSignificance < iJointSignificance) {
			optimalJointSignificance = iJointSignificance;
			optimalBoundaryBin = i;
		}
	}

	Double_t sIntegralBelow = _signal->Integral(0, optimalBoundaryBin);
	Double_t bIntegralBelow = _background->Integral(0, optimalBoundaryBin);
	Double_t sIntegralAbove = _signal->Integral(optimalBoundaryBin+1, _signal->GetNbinsX()+1);
	Double_t bIntegralAbove = _background->Integral(optimalBoundaryBin+1, _signal->GetNbinsX()+1);
	Double_t optimalSignificanceBelow = sIntegralBelow/std::sqrt(sIntegralBelow + bIntegralBelow);
	Double_t optimalSignificanceAbove = sIntegralAbove/std::sqrt(sIntegralAbove + bIntegralAbove);
	Double_t uncategorizedSignificance = _signal->Integral(0,_signal->GetNbinsX()+1)/std::sqrt(_signal->Integral(0,_signal->GetNbinsX()+1)+_background->Integral(0,_signal->GetNbinsX()+1));

	std::cout<<"**************************************************************************************************************************"<<std::endl<<
	"Optmal boundary for categorizing in the variable "<<_signal->GetXaxis()->GetTitle()<<std::endl<<
	"\t\tBin # = "<<optimalBoundaryBin<<std::endl<<
	"\t\tUp edge = "<<_signal->GetXaxis()->GetBinUpEdge(optimalBoundaryBin)<<std::endl<<
	"\t\tUncategorized Significance = "<<uncategorizedSignificance<<std::endl<<
	"\t\tOptimized Significance Above = "<<optimalSignificanceAbove<<std::endl<<
	"\t\tOptimized Significance Below = "<<optimalSignificanceBelow<<std::endl<<
	"\t\tJoint Significance = "<<std::sqrt(optimalSignificanceBelow*optimalSignificanceBelow + optimalSignificanceAbove*optimalSignificanceAbove)<<
	"**************************************************************************************************************************"<<std::endl;

	return _signal->GetXaxis()->GetBinUpEdge(optimalBoundaryBin);
};

std::vector<Double_t> vecString2vecDouble(std::vector<std::string> _numStrings) {
	std::vector<Double_t> _convertedDoubles;
	for (UInt_t i = 0; i < _numStrings.size(); i++) {
		std::string _iString =  _numStrings[i];
		trim(_iString);
		_convertedDoubles.push_back(std::stod(_iString));
	}
	return _convertedDoubles;
};


Bool_t stringIsInteger(std::string _isThisAnInt) {
	std::string::const_iterator it = _isThisAnInt.begin();
	while (it != _isThisAnInt.end() && std::isdigit(*it)) ++it;
	return !_isThisAnInt.empty() && it == _isThisAnInt.end();
};

Bool_t stringIsNumber(std::string _isThisANumber) {
	Bool_t is_a_number = false;

	try {
		boost::lexical_cast<double>(_isThisANumber);
		is_a_number = true;
	} catch (boost::bad_lexical_cast &) {}

	return is_a_number;
};


TDirectory *mkTFileDir(TFile *_file, std::string _dir) {
	if (_file->GetDirectory(_dir.c_str()) == nullptr) _file->mkdir(_dir.c_str());
	return _file->GetDirectory(_dir.c_str());
};


int caselessSubstringPos( std::string str1, std::string str2) {
	std::string::const_iterator it = std::search( str1.begin(), str1.end(),
		str2.begin(), str2.end(), [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2);});
	if ( it != str1.end() ) return it - str1.begin();
	else return -1; // not found
};


Bool_t findStringIC(const std::string & strHaystack, const std::string & strNeedle) {
	auto it = std::search(
		strHaystack.begin(), strHaystack.end(),
		strNeedle.begin(),   strNeedle.end(),
		[](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
		);
	return (it != strHaystack.end() );
};


Bool_t stringIsEmpty(std::string str) {
	if (str.find_first_not_of(' ') != std::string::npos)	return 0;
	return 1;
};


//line=style,hexcolor-alpha,width
template<typename T>
UChar_t setLineAtts(T* graph, std::string _atts, std::string _delimiter) {
	Int_t returnVal = 0;

	indexer<std::string> sequence({"style", "hexcolor-alpha", "width"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Line attributes of " + (std::string) graph->GetName() + " set to:";

	for (Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++) {
		if ((iAtt == sequence["style"]) && stringIsNumber(atts[iAtt])) {
			graph->SetLineStyle(std::stoi(atts[iAtt]));
			std::cout<<"\tstyle = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["hexcolor-alpha"]) && !stringIsEmpty(atts[iAtt])) {
			Int_t doTransparency = caselessSubstringPos( atts[iAtt], "-");
			if (doTransparency > -1 && stringIsNumber(atts[iAtt].substr(doTransparency+1))) {
				Float_t alpha = std::stof(atts[iAtt].substr(doTransparency+1));
				graph->SetLineColorAlpha(TColor::GetColor(atts[iAtt].c_str()), alpha);
				std::cout<<" color = " + atts[iAtt].substr(0,doTransparency-1) + " (alpha=" + atts[iAtt].substr(doTransparency+1) + ") ";
				setBit(returnVal,iAtt,1);
			} else {
				graph->SetLineColor(TColor::GetColor(atts[iAtt].c_str()));
				std::cout<<" color = " + atts[iAtt] + " ";
				setBit(returnVal,iAtt,1);
			}
		} else if ((iAtt == sequence["width"]) && stringIsNumber(atts[iAtt])) {
			graph->SetLineWidth(std::stof(atts[iAtt]));
			std::cout<<" width = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
	}

	std::cout<<std::endl;

	return returnVal;
};


//marker=style,hexcolor-alpha,size
template<typename T>
UChar_t setMarkerAtts(T* graph, std::string _atts, std::string _delimiter) {
	Int_t returnVal = 0;

	indexer<std::string> sequence({"style", "hexcolor-alpha", "size"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Marker attributes of " + (std::string) graph->GetName() + " set to:";

	for (Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++) {
		if ((iAtt == sequence["style"]) && stringIsNumber(atts[iAtt])) {
			graph->SetMarkerStyle(std::stoi(atts[iAtt]));
			std::cout<<"\tstyle = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["hexcolor-alpha"]) && !stringIsEmpty(atts[iAtt])) {
			Int_t doTransparency = caselessSubstringPos( atts[iAtt], "-");
			if (doTransparency > -1 && stringIsNumber(atts[iAtt].substr(doTransparency+1))) {
				Float_t alpha = std::stof(atts[iAtt].substr(doTransparency+1));
				graph->SetMarkerColorAlpha(TColor::GetColor(atts[iAtt].c_str()), alpha);
				std::cout<<" color = " + atts[iAtt].substr(0,doTransparency-1) + " (alpha=" + atts[iAtt].substr(doTransparency+1) + ") ";
				setBit(returnVal,iAtt,1);
			} else {
				graph->SetMarkerColor(TColor::GetColor(atts[iAtt].c_str()));
				std::cout<<" color = " + atts[iAtt] + " ";
				setBit(returnVal,iAtt,1);
			}
		} else if ((iAtt == sequence["width"]) && stringIsNumber(atts[iAtt])) {
			graph->SetMarkerSize(std::stof(atts[iAtt]));
			std::cout<<" size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
	}

	std::cout<<std::endl;

	return returnVal;
};


//fill=style,hexcolor-alpha
template<typename T>
UChar_t setFillAtts(T* graph, std::string _atts, std::string _delimiter) {
	Int_t returnVal = 0;

	indexer<std::string> sequence({"style", "hexcolor-alpha"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Fill attributes of " + (std::string) graph->GetName() + " set to:";


	for (Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++) {
		if ((iAtt == sequence["style"]) && stringIsNumber(atts[iAtt])) {
			graph->SetFillStyle(std::stoi(atts[iAtt]));
			std::cout<<"\tstyle = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["hexcolor-alpha"]) && !stringIsEmpty(atts[iAtt])) {
			Int_t doTransparency = caselessSubstringPos( atts[iAtt], "-");
			if (doTransparency > -1 && stringIsNumber(atts[iAtt].substr(doTransparency+1))) {
				Float_t alpha = std::stof(atts[iAtt].substr(doTransparency+1));
				graph->SetFillColorAlpha(TColor::GetColor(atts[iAtt].c_str()), alpha);
				std::cout<<" color = " + atts[iAtt].substr(0,doTransparency-1) + " (alpha=" + atts[iAtt].substr(doTransparency+1) + ") ";
				setBit(returnVal,iAtt,1);
			} else {
				graph->SetFillColor(TColor::GetColor(atts[iAtt].c_str()));
				std::cout<<" color = " + atts[iAtt] + " ";
				setBit(returnVal,iAtt,1);
			}
		}
	}

	std::cout<<std::endl;

	return returnVal;
};


//axis=range(min-max),limits(min-max),title_size,label_size,ndivs(n1-n2-n3-opt),title_offset,label_offset,more_log_labels,maxDigits
UChar_t setAxisAtts(TAxis* _axis, std::string _atts, std::string _delimiter) {
	Int_t returnVal = 0;

	indexer<std::string> sequence({"range(min-max)", "limits(min-max)", "title_size", "label_size", "ndivs", "title_offset", "label_offset", "more_log_labels", "maxDigits"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);

	std::cout<<"TAxis attributes of " + (std::string) _axis->GetName();
	if (_axis->GetParent()) std::cout<< " (" + (std::string) _axis->GetParent()->GetName() + ")";
	std::cout<< " set to:";

	for (Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++) {
		if ((iAtt == sequence["range(min-max)"]) && !stringIsEmpty(atts[iAtt])) {
			std::vector<Float_t> rangeAtts = strToFloatList(atts[iAtt], "-");
			if (rangeAtts.size() == 2) {
				_axis->SetRangeUser(rangeAtts[0], rangeAtts[1]);
				std::cout<<" range = " + atts[iAtt] + " ";
			}
		} else if ((iAtt == sequence["limits(min-max)"]) && !stringIsEmpty(atts[iAtt])) {
			std::vector<Float_t> limitsAtts = strToFloatList(atts[iAtt], "-");
			if (limitsAtts.size() == 2) {
				_axis->SetLimits(limitsAtts[0], limitsAtts[1]);
				std::cout<<" limits = " + atts[iAtt] + " ";
			}
		} else if ((iAtt == sequence["title_size"]) && stringIsNumber(atts[iAtt])) {
			_axis->SetTitleSize(std::stof(atts[iAtt]));
			std::cout<<" title size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["label_size"]) && stringIsNumber(atts[iAtt])) {
			_axis->SetLabelSize(std::stof(atts[iAtt]));
			std::cout<<" label size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["ndivs"]) && !stringIsEmpty(atts[iAtt])) {
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");
			std::vector<std::string> nDivAtt = split_string(atts[iAtt], "-");
			_axis->SetNdivisions(std::stoi(nDivAtt[0]), std::stoi(nDivAtt[1]), std::stoi(nDivAtt[2]), (nDivAtt.size()>3) ? nDivAtt[3].c_str() : "");
			std::cout<<" nDivisions = " + nDivAtt[0] + " " + nDivAtt[1] + " " + nDivAtt[2] + " " + ((nDivAtt.size()>3) ? nDivAtt[3] : "") + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["title_offset"]) && stringIsNumber(atts[iAtt])) {
			_axis->SetTitleOffset(std::stof(atts[iAtt]));
			std::cout<<" title offset = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["label_offset"]) && stringIsNumber(atts[iAtt])) {
			_axis->SetLabelOffset(std::stof(atts[iAtt]));
			std::cout<<" label offset = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["more_log_labels"]) && stringIsNumber(atts[iAtt])) {
			_axis->SetMoreLogLabels(std::stoi(atts[iAtt]));
			std::cout<<" more log labels = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["maxDigits"]) && stringIsNumber(atts[iAtt])) {
			((TGaxis*)_axis)->SetMaxDigits(std::stoi(atts[iAtt]));
			std::cout<<" max digits = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}

	}

	std::cout<<std::endl;

	return returnVal;
};


// pad=coordinates(x1-y1-x2-y2),margins(L-R-B-T),fillStyle-color-alpha,bordersize,grid(x-y),log(x-y)
template<typename T>
UChar_t setPadAtts(T* _pad, std::string _atts, std::string _delimiter) {
	Int_t returnVal = 0;

	indexer<std::string> sequence({"coordinates(x1-y1-x2-y2)", "margins(L-R-B-T)", "fillStyle-color-alpha", "bordersize", "grid(x-y)", "log(x-y)"});

	_atts = findAndReplaceAll(_atts, " ", "");
	std::vector<std::string> atts = split_string(_atts, _delimiter);
	std::cout<<"Pad attributes of " + (std::string) _pad->GetName() + " set to:";

	for (Int_t iAtt = 0; (iAtt < (Int_t) atts.size()) && (iAtt < (Int_t) sequence.size()); iAtt++) {
		if ((iAtt == sequence["coordinates(x1-y1-x2-y2)"]) && !stringIsEmpty(atts[iAtt])) {
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");
			std::vector<std::string> coords = split_string(atts[iAtt], "-");
			if (coords.size() == 4 && stringIsNumber(coords[0]) && stringIsNumber(coords[1]) && stringIsNumber(coords[2]) && stringIsNumber(coords[3])) {
				_pad->SetPad(std::stof(coords[0]), std::stof(coords[1]), std::stof(coords[2]), std::stof(coords[3]));
				std::cout<<" (x1,y1,x2,y2) = (" + coords[0] + "," + coords[1] + "," + coords[2] + "," + coords[2] + ")";
				setBit(returnVal,iAtt,1);
			}
		} else if ((iAtt == sequence["margins(L-R-B-T)"]) && !stringIsEmpty(atts[iAtt])) {
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");
			std::vector<std::string> margins = split_string(atts[iAtt], "-");
			if (stringIsNumber(margins[0]) && stringIsNumber(margins[1]) && stringIsNumber(margins[2]) && stringIsNumber(margins[3])) {
				_pad->SetMargin(std::stof(margins[0]), std::stof(margins[1]), std::stof(margins[2]), std::stof(margins[3]));
				std::cout<<" margins(L,R,B,T) = (" + margins[0] + "," + margins[1] + "," + margins[2] + "," + margins[3] + ")";
				setBit(returnVal,iAtt,1);
			}
		} else if ((iAtt == sequence["fillStyle-color-alpha"]) && !stringIsEmpty(atts[iAtt])) {
			setFillAtts(_pad, atts[iAtt], "-");
		} else if ((iAtt == sequence["bordersize"]) && stringIsNumber(atts[iAtt])) {
			_pad->SetBorderSize(std::stof(atts[iAtt]));
			std::cout<<" border size = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		} else if ((iAtt == sequence["grid(x-y)"]) && !stringIsEmpty(atts[iAtt])) {
			atts[iAtt] = findAndReplaceAll(atts[iAtt], " ", "");

			std::vector<std::string> gridAtts = split_string(atts[iAtt], "-");

			if (gridAtts.size() > 0 && stringIsNumber(gridAtts[0])) {
				_pad->SetGridx(std::stoi(gridAtts[0]));
				std::cout<<" grid X = " + gridAtts[0] + " ";
				setBit(returnVal,iAtt,1);
			}

			if (gridAtts.size() > 1 && stringIsNumber(gridAtts[1])) {
				_pad->SetGridy(std::stoi(gridAtts[1]));
				std::cout<<" grid Y = " + gridAtts[1] + " ";
				setBit(returnVal,iAtt,1);
			}
		} else if ((iAtt == sequence["log(x-y)"]) && !stringIsEmpty(atts[iAtt])) {
			std::vector<std::string> logXY = split_string(atts[iAtt], "-");
			if (logXY.size()>0 && stringIsNumber(logXY[0])) _pad->SetLogx(std::stoi(logXY[0]));
			if (logXY.size()>1 && stringIsNumber(logXY[1])) _pad->SetLogy(std::stoi(logXY[1]));
			std::cout<<" log scale (X,Y) = " + atts[iAtt] + " ";
			setBit(returnVal,iAtt,1);
		}
	};

	std::cout<<std::endl;

	return returnVal;
};


// frame=fill=style,color-alpha;line=style,color-alpha,width
template<typename T>
void setFrameAtts(T* _hist, std::string _style, std::string _del1, std::string _del2) {
	_style = findAndReplaceAll(_style, " ", "");
	std::vector<std::string> styleOptions = split_string(_style, _del1);
	for (const std::string styleOption : styleOptions) {
		cout<<styleOption<<endl;
		std::vector<std::string> varVal = split_string(styleOption, "=");
		if (stringToUpper(varVal[0]) == stringToUpper("line")) setLineAtts(_hist, varVal[1], _del2);
		else if (stringToUpper(varVal[0]) == stringToUpper("fill")) setFillAtts(_hist, varVal[1], _del2);
	};
};


//histstyle==line=style,hexcolor-alpha,width;marker=style,hexcolor-alpha,size;fill=style,hexcolor-alpha
template<typename T>
void setHistStyle(T* _hist, std::string _style, std::string _del1, std::string _del2) {
	_style = findAndReplaceAll(_style, " ", "");
	std::vector<std::string> styleOptions = split_string(_style, _del1);
	for (const std::string styleOption : styleOptions) {
		std::vector<std::string> varVal = split_string(styleOption, "=");
		if (stringToUpper(varVal[0]) == stringToUpper("marker")) setMarkerAtts(_hist, varVal[1], _del2);
		else if (stringToUpper(varVal[0]) == stringToUpper("line")) setLineAtts(_hist, varVal[1], _del2);
		else if (stringToUpper(varVal[0]) == stringToUpper("fill")) setFillAtts(_hist, varVal[1], _del2);
	};
};


// legendstyle==bounds=border_size,x1-y1-x2-y2;fill=style,hexcolor-alpha;entries=text_size,nCols
void setLegendStyle(TLegend *_legend, std::string _style, std::string _del1=";", std::string _del2=",") {
	std::vector<std::string> legendStyleOptions = split_string(_style, _del1);

	std::cout<<"Legend attributes of "<<_legend->GetName()<< " set to ";
	for (const std::string styleOption : legendStyleOptions) {
		std::vector<std::string> subOptions = split_string(styleOption, "=");

		if (stringToUpper(subOptions[0]) == stringToUpper("bounds")) {
			std::vector<std::string> boundsOpts = split_string(subOptions[1], _del2);
			if (stringIsNumber(boundsOpts[0])) {
				_legend->SetBorderSize(std::stof(boundsOpts[0]));
				std::cout<<"border size = "+boundsOpts[0]<<" ";
			}
			std::vector<Float_t> coordinates = (boundsOpts.size()>0) ? strToFloatList(boundsOpts[1],"-") : (std::vector<Float_t>) {};
			if (coordinates.size() == 4) {
				_legend->SetX1(coordinates[0]);
				_legend->SetY1(coordinates[1]);
				_legend->SetX2(coordinates[2]);
				_legend->SetY2(coordinates[3]);
				std::cout<<"coordinates = "+boundsOpts[1]<<" ";
			}
		} else if (stringToUpper(subOptions[0]) == stringToUpper("fill")) {
			setFillAtts(_legend, subOptions[1], _del2);
		} else if (stringToUpper(subOptions[0]) == stringToUpper("entries")) {
			std::vector<std::string> entryOpts = split_string(subOptions[1], _del2);
			if (entryOpts.size()>0 && stringIsNumber(entryOpts[0])) {
				_legend->SetTextSize(std::stof(entryOpts[0]));
				std::cout<<"text size = "+entryOpts[0]<<" ";
			}
			if (entryOpts.size()>1 && stringIsNumber(entryOpts[1])) {
				_legend->SetNColumns(std::stoi(entryOpts[1]));
				std::cout<<"N columns = "+entryOpts[1]<<" ";
			}
		}
	};
	std::cout<<std::endl;
};


struct DinkyHister {

	DinkyHister() {};

	~DinkyHister() {
		for (TH1 * hist : histograms) {
			hist->Delete();
		};
	};

	DinkyHister(std::string _canvasSize="1600x1200", std::string _padAtts="", std::string _legendStyle = "bounds=0.,0.8-0.7-0.97-0.95;fill=4000,;entries=0.1,1",
		std::string _titles="") {
		std::vector<std::string> canvasSize = split_string(_canvasSize,"x");
		canvas.SetCanvasSize(1500, 1500);
		canvas.SetWindowSize(std::stoi(canvasSize[0]), std::stoi(canvasSize[1]));
		setPadAtts(&canvas, _padAtts);
		setLegendStyle(&legend, _legendStyle);
		histStack.SetTitle(_titles.c_str());
	};

	void applyStyleToAll(std::string _style, std::string _del1=";", std::string _del2=",") {
		for (TH1* hist : histograms) {
			setHistStyle(hist, _style, _del1, _del2);
		}
	};
	//_legend=legend--option
	void add(TH1* _hist, std::string _style, std::string _drawOption, std::string _legend) {

		std::string newName = (std::string) _hist->GetName() + "_" + std::to_string(reinterpret_cast<ULong_t>(this)) + "_" + std::to_string(reinterpret_cast<ULong_t>(_hist));

		TH1* hist = (TH1*) _hist->Clone(newName.c_str());
		setHistStyle(hist, _style);
		histStack.Add(hist, _drawOption.c_str());
		histograms.push_back(hist);

		if (!stringIsEmpty(_legend)) {
			std::vector<std::string> legendOpt = split_string(_legend, "--");
			TLegendEntry* legendEntry = nullptr;
			if (legendOpt.size()>1) legendEntry = legend.AddEntry(hist, legendOpt[0].c_str(), legendOpt[1].c_str());
			else legendEntry = legend.AddEntry(_hist, legendOpt[0].c_str());
			legendEntry->SetTextColor(hist->GetLineColor());
		}
	};

	void draw(std::string _option) {
		canvas.Draw();
		canvas.cd();
		histStack.Draw(_option.c_str());
		legend.Draw();
	};

	void write(std::string _path, std::string _fileName, std::string _formats = "png,pdf") {
		_formats = findAndReplaceAll(_formats, ".", "");
		std::vector<std::string> extensions = split_string(_formats, ",");
		mkdir(_path);
		for (const std::string & iFormat : extensions) {
			if (!stringIsEmpty(iFormat)) canvas.SaveAs((_path + "/" + _fileName + "." + iFormat).c_str());
		}

	};

	void disown() {
		histograms.clear();
	};

	void update() {
		gPad->Update();
		gPad->RedrawAxis();
		gPad->Modified();

		canvas.Update();
		canvas.RedrawAxis();
		canvas.Modified();
	};

	std::vector<TH1*> histograms;
	TCanvas canvas;
	TLegend legend;
	THStack histStack;
	std::string titles;
};


struct correlationMatix {
	UInt_t 								nVar;
	std::vector<std::vector<Double_t>>	data2p;
	std::vector<Double_t>				data1p;
	Double_t 							sumW;
	Bool_t 								doMagic = 0;
	Double_t 							magicFactor;
	Double_t 							magicFactor2;
	UInt_t  							nanCounter = 0;

	correlationMatix() {};

	correlationMatix(Int_t _nVar, Bool_t _doMagic = 0, Double_t _magicFactor = 1000) {
		nVar 			=	_nVar;
		data2p			=	std::vector<std::vector<Double_t>>(nVar, std::vector<Double_t>(nVar, 0.));
		data1p			=	std::vector<Double_t>(nVar, 0.);
		sumW 			= 	0.;
		doMagic 		= 	_doMagic;
		magicFactor     = _magicFactor;
		magicFactor2 	= magicFactor * magicFactor;
	};

	void addPoint(std::vector<Double_t> _point, Double_t _weight) {
		if (std::isinf( _weight) || std::isnan(_weight)) return;

		UInt_t summedCounter = 0;

		for (UInt_t i = 0; i < nVar; i++) {

			if (std::isinf(_point[i]) || std::isnan(_point[i])) {
				nanCounter += 1;
				continue;
			}

			UInt_t counteri = 0;

			for (UInt_t j = 0; j < nVar; j++) {

				if (std::isinf(_point[j]) || std::isnan(_point[j])) {
					nanCounter += 1;
					continue;
				}

				if (doMagic) data2p[i][j] += (_point[i]*_point[j]*_weight) / magicFactor2;
				else data2p[i][j] += (_point[i]*_point[j]*_weight);

				counteri += 1;
			}

			if (counteri > 0) {

				if (doMagic) data1p[i] += _point[i]*_weight / magicFactor;
				else data1p[i] += _point[i]*_weight;

				summedCounter += 1;
			}
		}

		if (summedCounter > 0) sumW += _weight;
	}

	TH2F* getCorrelationHistAbs() {
		TH2F* 						corrHist 		=		new TH2F("corrHist", "", nVar, 0, nVar, nVar, 0, nVar);
		std::vector<Double_t>		means;

		for (UInt_t i = 0; i < data1p.size(); i++) {
			means.push_back(data1p[i]/sumW);
		}

		for (UInt_t i = 0; i < data2p.size(); i++) {
			Double_t 	denI 	= std::sqrt(data2p[i][i] - sumW * means[i]* means[i]);

			for (UInt_t j = 0; j < data2p.size(); j++) {

				Double_t  Rij	=	data2p[i][j] - sumW * means[i]*means[j];
				Rij 			/=	denI;
				Rij 			/=	std::sqrt(data2p[j][j] - sumW * means[j]* means[j]);

				corrHist->SetBinContent(1+i, 1+j, std::abs(Rij));
			}
		}

		return corrHist;
	}

	TH2F* getCorrelationHist() {
		TH2F* 						corrHist 		=		new TH2F("corrHist", "", nVar, 0, nVar, nVar, 0, nVar);
		std::vector<Double_t>		means;

		for (UInt_t i = 0; i < data1p.size(); i++) {
			means.push_back(data1p[i]/sumW);
		}

		for (UInt_t i = 0; i < data2p.size(); i++) {
			Double_t 	denI 	= std::sqrt(data2p[i][i] - sumW * means[i]* means[i]);
			for (UInt_t j = 0; j < data2p.size(); j++) {
				Double_t  Rij	=	data2p[i][j] - sumW * means[i]*means[j];
				Rij 			/=	denI;
				Rij 			/=	std::sqrt(data2p[j][j] - sumW * means[j]* means[j]);
				std::cout<<"r["<<i<<", "<<j<<"]\t=\t"<<Rij<<std::endl;
				corrHist->SetBinContent(1+i, 1+j, Rij);
			}
		}

		std::cout<<"NAN occurrences: "<<nanCounter<<std::endl;

		return corrHist;
	}

};


std::string stringToUpper(std::string _str) {
	boost::to_upper(_str);
	return _str;
};


std::string stringToLower(std::string _str) {
	boost::to_lower(_str);
	return _str;
};


template<typename T1, typename T2>
Int_t findIndex(const std::vector<T1> & _haystack, T2 _needle) {
	Int_t Index = -999;
	for (UInt_t i = 0; i < _haystack.size(); i++) {
		if ( _haystack[i] == (T1) _needle) {
			Index = i;
			break;
		}
	}
	return Index;
};


Bool_t objectExists(TFile *_File, std::string _objectName) {
	TObject *_object = _File->Get(_objectName.c_str());
	if (!_object) return 0;
	_object->Delete();
	return 1;
};


Bool_t objectExists(std::string _FilePath, std::string _objectName) {
	if (!file_exists(_FilePath)) return 0;
	TFile _file(_FilePath.c_str(), "READ");
	TObject *_object = _file.Get(_objectName.c_str());
	if (!_object) return 0;
	_object->Delete();
	_file.Close();
	return 1;
};


std::vector<Float_t> strToFloatList(std::string _listString, std::string _delimiter) {
	std::vector<std::string> numberStrings = split_string(_listString, _delimiter);
	std::vector<Float_t> convertedNumbers;
	for (const std::string & numberString : numberStrings) {
		if (stringIsEmpty(numberString)) convertedNumbers.push_back(std::numeric_limits<Float_t>::max());
		else if (stringIsNumber(numberString)) convertedNumbers.push_back(std::stof(numberString));
	};
	return convertedNumbers;
};

std::vector<Double_t> strToDoubleList(std::string _listString, std::string _delimiter) {
	std::vector<std::string> numberStrings = split_string(_listString, _delimiter);
	std::vector<Double_t> convertedNumbers;
	for (const std::string & numberString : numberStrings) {
		if (stringIsEmpty(numberString)) convertedNumbers.push_back(std::numeric_limits<Double_t>::max());
		else if (stringIsNumber(numberString)) convertedNumbers.push_back(std::stod(numberString));
	};
	return convertedNumbers;
};


Int_t hex2rootColor(std::string _hexCode) {
	trim(_hexCode);
	if (!match("#*",_hexCode)) _hexCode = "#" + _hexCode;
	return TColor::GetColor(_hexCode.c_str());
};


std::string sysexec(std::string cmd) {
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
	if (!pipe) {
		throw std::runtime_error("popen() failed!");
	}
	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
		result += buffer.data();
	}
	return result;
};


std::vector<Double_t> getNpercentMinInterval(TH1 *_hist, Float_t _N) {
	if (_hist->IsZombie()) {
		cout << "Error: Null histogram!" << endl;
		exit (EXIT_FAILURE);
	}

	UInt_t Nbins = _hist->GetNbinsX();
	UInt_t firstBin = 1, lastBin = Nbins;
	Float_t P = _N / 100.0;
	TH1D* cum = (TH1D*) _hist->GetCumulative(1, "cum");
	TH1D* rcum = (TH1D*) _hist->GetCumulative(0, "rcum");
	TAxis *xaxis = _hist->GetXaxis();
	Double_t integral = cum->GetBinContent(Nbins);
	Double_t tempInterval = (xaxis->GetXmax()) - (xaxis->GetXmin());
	Double_t partialIntegral = 0;
	Double_t leftEdge, rightEdge;

	for (UInt_t i = 1; i <= Nbins; i++) {
		if (P * integral > rcum->GetBinContent(i)) break;
		leftEdge = xaxis->GetBinLowEdge(i);
		for (UInt_t j = i; j < Nbins + 1; j++) {
			partialIntegral = cum->GetBinContent(j) - cum->GetBinContent(i - 1);
			rightEdge = xaxis->GetBinUpEdge(j);
			if ((partialIntegral / integral) >= P && (tempInterval > (rightEdge - leftEdge))) {
				lastBin = j;
				firstBin = i;
				tempInterval = rightEdge - leftEdge;
				break;
			}
		}
	}
	cum->Delete();
	rcum->Delete();

	std::vector<Double_t> minRange = {xaxis->GetBinLowEdge(firstBin), xaxis->GetBinUpEdge(lastBin)};

	return minRange;
};



struct logStream {

	logStream() {};

	logStream(std::string _logFilePath) {
		init(_logFilePath);
	};

	~logStream() {
		logFile.close();
	};

	void init(std::string _logFilePath) {
		logFile.open(_logFilePath);
		std::cout<< "Writing logs to file: "<<_logFilePath<<std::endl;
	};

	ofstream logFile;

	template<typename T>
	logStream & operator<<(const T& mValue) {
		std::cout << mValue;
		logFile << mValue;
		return *this;
	};
};


// std::vector<Double_t> getNpercentLinearInterval(TH1 *_hist, Float_t _N){
// 	return {};
// };


std::vector<std::string> prefixVecString(std::vector<std::string> _vecStr, std::string _prefix, Int_t _insPos) {
	std::vector<std::string> prefixedVecStr = _vecStr;
	if (_insPos < 0) {
		std::cout<<"Error! negative insert position passed to prefixVecString()"<<std::endl;
		return {};
	}

	UInt_t tmpInsPos = _insPos;
	if (_insPos < 0) {
		std::for_each(prefixedVecStr.begin(), prefixedVecStr.end(), [&](std::string & iString) {
			iString.append(_prefix);
		});
	} else {
		std::for_each(prefixedVecStr.begin(), prefixedVecStr.end(), [&](std::string & iString) {
			if (_insPos > (Int_t) iString.length())  tmpInsPos = iString.length();
			iString.insert(tmpInsPos, _prefix);
		});
	}
	return prefixedVecStr;
};



void normalizeHist(TH1* _hist, Double_t _norm, Bool_t _full) {
	if (_full)_hist->Scale(_norm/_hist->Integral(0, _hist->GetNbinsX()+1));
	else _hist->Scale(_norm/_hist->Integral());
};


void normalizeHist(TH1& _hist, Double_t _norm, Bool_t _full) {
	if (_full)_hist.Scale(_norm/_hist.Integral(0, _hist.GetNbinsX()+1));
	else _hist.Scale(_norm/_hist.Integral());
};


std::vector<Int_t> strToIntList(std::string _listString, std::string _delimiter) {
	std::vector<std::string> numberStrings = split_string(_listString, _delimiter);
	std::vector<Int_t> convertedNumbers;
	for (const std::string & numberString : numberStrings) {
		if (stringIsEmpty(numberString)) convertedNumbers.push_back(std::numeric_limits<Int_t>::max());
		else if (stringIsNumber(numberString)) convertedNumbers.push_back(std::stoi(numberString));
	};
	return convertedNumbers;
};


Double_t spearmanR(TH2 *_hist, Int_t _Nbins) {

	TH1D*_xProj = (TH1D*)_hist->ProjectionX();
	normalizeHist(_xProj, 1.);
	TH1D*_xCum = (TH1D*)_xProj->GetCumulative();
	delete _xProj;

	TH1D*_yProj = (TH1D*)_hist->ProjectionY();
	normalizeHist(_yProj, 1.);
	TH1D*_yCum = (TH1D*)_yProj->GetCumulative();
	delete _yProj;

	std::string xyRanksName = (std::string)_hist->GetName() + "_xyRanks";
	TH2D xyRanks(xyRanksName.c_str(), "", _Nbins, 0., 1.000001, _Nbins,  0., 1.000001);

	Int_t _xNbins = _xCum->GetNbinsX() + 1;
	Int_t _yNbins = _yCum->GetNbinsX() + 1;

	for (Int_t iX = 1; iX < _xNbins; iX++) {
		Double_t iXrank = _xCum->GetBinContent(iX);
		for (Int_t iY = 1; iY < _yNbins; iY++) {
			Double_t iYrank = _yCum->GetBinContent(iY);
			Double_t iXYweight = _hist->GetBinContent(iX, iY);
			xyRanks.Fill(iXrank, iYrank, iXYweight);
		}
	}

	delete _xCum;
	delete _yCum;

	Double_t spearmanRval = xyRanks.GetCorrelationFactor();

	return spearmanRval;
};


Float_t relInvMass(Float_t _pT1, Float_t _eta1, Float_t _phi1, Float_t _pT2, Float_t _eta2, Float_t _phi2) {
	////	https://en.wikipedia.org/wiki/Invariant_mass#Collider_experiments
	Float_t invMassVal = std::sqrt(2.*_pT1*_pT2*(std::cosh(_eta1-_eta2) - std::cos(_phi1-_phi2)));
	return invMassVal;
};



Float_t relTransMass(Float_t _pT1, Float_t _phi1, Float_t _pT2, Float_t _phi2) {
	if (_pT1 < 0. || _pT2 < 0.) {
		std::cout<<"Error! Negative pT given!!!!!!!!!!"<<std::endl;
		return -9999.;
	}
	Float_t relTransMassVal 	=	std::sqrt(2. * _pT1 * _pT2 * (1. - std::cos(_phi1 - _phi2)));
	return relTransMassVal;
};



template<typename T>
void remove_intersection(std::vector<T>& a, std::vector<T>& b) {
	std::unordered_multiset<T> st;
	st.insert(a.begin(), a.end());
	st.insert(b.begin(), b.end());
	auto predicate = [&st](const T& k) { return st.count(k) > 1; };
	a.erase(std::remove_if(a.begin(), a.end(), predicate), a.end());
	b.erase(std::remove_if(b.begin(), b.end(), predicate), b.end());
};


TH1* getSqrtHist(TH1* _hist) {
	std::string sqrtName = _hist->GetName();
	sqrtName += "_sqrt";
	TH1* _sqrtHist = (TH1*) _hist->Clone(sqrtName.c_str());
	_sqrtHist->Reset();

	for (Int_t iBinX = 1; iBinX <= _hist->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= _hist->GetNbinsY(); iBinY++) {

			Double_t 					binContent		=		_hist->GetBinContent(iBinX, iBinY);
			binContent 									=		(binContent>0) ? std::sqrt(binContent) : 0.;								///// 0 if negative
			Double_t 					binErr 			=		(binContent>0) ? 0.5 * _hist->GetBinError(iBinX, iBinY)/binContent : 0.;	//// sigma(x^0.5) = 0.5 * sigma(x)/x^0.5
			_sqrtHist->SetBinContent(iBinX, iBinY, binContent);
			_sqrtHist->SetBinError(iBinX, iBinY, binErr);
		}
	}

	return _sqrtHist;
};


Bool_t hasNegativeBins(TH1* _hist, Float_t _thres) {

	Bool_t tmpHasNegativeBins = 0;
	std::string tmpStr = std::string("Negative Bins found in hist ") + std::string(_hist->GetName()) + " :\n";

	for (Int_t iBinX = 1; iBinX <= _hist->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= _hist->GetNbinsX(); iBinY++) {
			Double_t iBinContent = _hist->GetBinContent(iBinX, iBinY);
			if (iBinContent < _thres) {
				tmpHasNegativeBins = 1;
				tmpStr += std::string("\t\t\t\t\t(") + std::to_string(iBinX) + std::string(", ") + std::to_string(iBinY) + std::string(") = ") +std::to_string(iBinContent) + "\n";
			}
		}
	}
	if (tmpHasNegativeBins) std::cout<<tmpStr;
	return tmpHasNegativeBins;
};


void removeNegativeBins(TH1* _hist, Float_t _thres) {
	for (Int_t iBinX = 1; iBinX <= _hist->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= _hist->GetNbinsX(); iBinY++) {
			Double_t iBinError = _hist->GetBinError(iBinX, iBinY, 0.);
			if (_hist->GetBinContent(iBinX, iBinY) < _thres) {
				_hist->SetBinContent(iBinX, iBinY, 0.);
				_hist->SetBinError(iBinX, iBinY, iBinError);
			}
		}
	}
};


TH1* copyHistSubRange(TH1* _hist, Float_t _xMin, Float_t _xMax, std::string newName) {

	Int_t xMinBin = _hist->GetXaxis()->FindBin(_xMin);
	Int_t xMaxBin = _hist->GetXaxis()->FindBin(_xMax);

	std::string copyName = newName;
	if(copyName.empty()) copyName = std::string(_hist->GetName()) + "_xRange_" + removeTrailingZeros(_xMin) + "_" + removeTrailingZeros(_xMax);
	copyName = findAndReplaceAll(copyName, ".", "p");

	TH1* hCopy = (TH1*) _hist->Clone(copyName.c_str());

	for (Int_t iBin = 1; iBin <= _hist->GetNbinsX(); iBin++) {
		if (iBin < xMinBin || iBin > xMaxBin) {
			hCopy->SetBinContent(iBin, 0.);
			hCopy->SetBinError(iBin, 0.);
		}
	}

	return hCopy;
};


bool isInteger(const std::string & s) {
	if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;
	char * p;
	strtol(s.c_str(), &p, 10);
	return (*p == 0);
};


Float_t foldedPhi(Float_t _phi) {
	if (std::abs(_phi) <= TMath::PiOver2()) return _phi;
	else if (_phi < - TMath::PiOver2()) return -(_phi + TMath::Pi());
	else if (_phi > TMath::PiOver2()) return -(_phi - TMath::Pi());
	return _phi;
};


Char_t RGB2ColPlt(std::string _CSVFile, Int_t _firstCol=0) {
	CSVReader csvFile(_CSVFile, ",");
	std::vector<std::vector<std::string>> csvDat = csvFile.getData();

	std::vector<Int_t> palette;

	for (UInt_t iRow = 0; iRow < csvDat.size(); iRow++) {
		if (csvDat[iRow].size() < UInt_t(_firstCol+3)) continue;
		if (!stringIsNumber(csvDat[iRow][_firstCol+0]) || !stringIsNumber(csvDat[iRow][_firstCol+1]) || !stringIsNumber(csvDat[iRow][_firstCol+2])) continue;
		if (isInteger(csvDat[iRow][_firstCol+0]) && isInteger(csvDat[iRow][_firstCol+1]) && isInteger(csvDat[iRow][_firstCol+2])) palette.push_back(TColor::GetColor(std::stoi(csvDat[iRow][_firstCol+0]), std::stoi(csvDat[iRow][_firstCol+1]), std::stoi(csvDat[iRow][_firstCol+2])));
		else palette.push_back(TColor::GetColor(std::stof(csvDat[iRow][_firstCol+0]), std::stof(csvDat[iRow][_firstCol+1]), std::stof(csvDat[iRow][_firstCol+2])));
	}

	gStyle->SetPalette(palette.size(), palette.data());

	std::cout<<"Loaded color palette from file "<<_CSVFile<<"\t\tN = "<<palette.size()<<std::endl;

	return 0;
};


TH1* getAbsHist(TH1* _hist) {
	TH1* _clone = (TH1*) _hist->Clone("abs_clone");
	_clone->Reset();

	for (Int_t iBinX = 1; iBinX <= _hist->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= _hist->GetNbinsY(); iBinY++) {
			_clone->SetBinContent(iBinX, iBinY, std::abs(_hist->GetBinContent(iBinX, iBinY)));
		}
	}
	return _clone;
};


template <typename T>
std::string ToSciString(T value, UShort_t nDec, std::string powFormat) {
	std::stringstream out;
	Int_t pow = std::floor(std::log10(std::abs(value)));
	Double_t tmpValue = value/std::pow(10,pow);
	if (pow == 0) {
		// if (std::abs(pow) >= nDec) {
		out << std::setprecision(nDec)<< std::fixed << tmpValue;
	} else {
		if (powFormat == "e") out << std::setprecision(nDec)<< std::scientific << value;
		else if (powFormat=="^" && pow==1)  out << std::setprecision(nDec)<<std::fixed<<tmpValue<<"#times10";
		else if (powFormat=="^") out << std::setprecision(nDec)<<std::fixed<<tmpValue<<"#times10^{"<<pow<<"}";
	}

	return out.str();
}


template <typename T>
std::string ToSciString(T value, T valErr, UShort_t nDec, std::string pmSign, std::string powFormat) {

	valErr = std::abs(valErr);
	Int_t minPow = std::min(std::floor(std::log10(std::abs(value))), std::floor(std::log10(valErr)));

	std::stringstream out;
	if (minPow==0 || minPow > -3) {
		out << std::setprecision(std::max(Int_t(nDec), std::abs(minPow)))<<std::fixed<<"("<<value<<pmSign<<valErr<<")";
	} else {

		value = value/std::pow(10,minPow);
		valErr = valErr/std::pow(10,minPow);

		if (powFormat == "e") out << std::setprecision(nDec)<<std::fixed<<"("<<value<<pmSign<<valErr<<")e"<<minPow;
		else if (powFormat=="^" && minPow==1) out << std::setprecision(nDec)<<std::fixed<<"("<<value<<pmSign<<valErr<<")#times10";
		else if (powFormat=="^") out << std::setprecision(nDec)<<std::fixed<<"("<<value<<pmSign<<valErr<<")#times10^{"<<minPow<<"}";
	}

	return out.str();
}


void cmsLegend(TPad* _pad, std::string _runInfo, std::string _2ndary, Float_t _scale) {

	/////// https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines
	///////	https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines/CMS_lumi.h
	/////// https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines/CMS_lumi.C
	static constexpr Float_t lumiOverCmsTextSize     	= 	0.8;
	static constexpr Float_t extraOverCmsTextSize  		= 	0.76;

	_pad->cd();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextFont(42); 	//// 4=helvetica-medium-r-normal, 2=scalable and rotatable hardware fonts

	Float_t ch = _pad->GetAbsHNDC();
	Float_t ctm = _pad->GetTopMargin();
	Float_t cmsTextSize = 0.98*ctm/ch;
	Float_t yPos = 1. - 0.75*ctm;

	latex.SetTextSize(cmsTextSize*_scale);
	latex.SetTextAlign(11);	//// bottom-left align
	latex.DrawLatex(_pad->GetLeftMargin(), yPos, "#bf{CMS} ");

	if (!_2ndary.empty()) {
		_2ndary = std::string("#it{  ") + _2ndary + "}";
		latex.SetTextSize(extraOverCmsTextSize*cmsTextSize*_scale);
		latex.DrawLatex(_pad->GetLeftMargin() + _scale*1.9*ctm, yPos, _2ndary.c_str());
	}

	if (!_runInfo.empty()) {
		latex.SetTextAlign(31);	//// bottom-right align
		latex.SetTextSize(lumiOverCmsTextSize*cmsTextSize*_scale);
		latex.DrawLatex(1.-_pad->GetRightMargin(), yPos, _runInfo.c_str());
	}
};


void TGraphErrorsToCSV(TGraphErrors *_graph, std::string _CSVoutfile) {

	std::string outStr="x, xErr, y, yErr\n";

	for (Int_t iPoint = 0; iPoint < _graph->GetN(); iPoint++) {
		Double_t iX, iY;
		_graph->GetPoint(iPoint, iX, iY);
		outStr += ToSciString(iX, 6) + ", " + ToSciString(_graph->GetErrorX(iPoint), 6) + ", " + ToSciString(iY, 6) + ", " + ToSciString(_graph->GetErrorY(iPoint), 6) + "\n";
	}

	ofstream frOutFile(_CSVoutfile);
	frOutFile << outStr;
	frOutFile.close();
	if (file_exists(_CSVoutfile)) std::cout<<"CSV written to "<<_CSVoutfile<<std::endl;
	else std::cout<<"Failed to write CSV to "<<_CSVoutfile<<std::endl;

};


void clearTHStack(THStack *hStack) {
	TIter 														hNext(hStack->GetHists());
	TObject*													iHist;
	while ((iHist = (TObject*) hNext())) {
		delete iHist;
	}
};


Float_t getPunziSmin(Float_t _B, Float_t _a, Float_t _b) {
	_B = (_B<0) ? 0. : _B;
	Float_t punziSminVal = _a*_a/8. + 	9.*_b*_b/13. + _a*std::sqrt(_B) + _b/2. * std::sqrt(_b*_b + 4.*_a*std::sqrt(_B) + 4.*_B);
	return punziSminVal;
};


TH1* getPunziSminHist(TH1* _Bhist, Float_t _a, Float_t _b) {

	TH1* hPunziSmin 							=			(TH1*) _Bhist->Clone((std::string(_Bhist->GetName())+"_Punzi").c_str());
	hPunziSmin->Reset();

	for (Int_t iBinX = 1; iBinX <= hPunziSmin->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= hPunziSmin->GetNbinsY(); iBinY++) {
			hPunziSmin->SetBinContent(iBinX,iBinY, getPunziSmin(_Bhist->GetBinContent(iBinX, iBinY), _a, _b));
		}
	}

	return hPunziSmin;
};


TGraphErrors* makeROC(TH1* _sigEffHist, TH1* _bgEffHist) {

	if ((_sigEffHist->GetNbinsX() != _bgEffHist->GetNbinsX()) || (_sigEffHist->GetNbinsY() != _bgEffHist->GetNbinsY())) {
		std::cout<<"Error making ROC! Unequal nbins!"<<std::endl;
		return nullptr;
	}

	TGraphErrors* 										sbROC 		= 		new TGraphErrors();
	sbROC->SetName((std::string(_sigEffHist->GetName()) + "_vs_" + std::string(_bgEffHist->GetName()) + "_ROC").c_str());

	for (Int_t iBinX = 1; iBinX <= _bgEffHist->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= _bgEffHist->GetNbinsX(); iBinY++) {
			addPointToGraph(sbROC, _bgEffHist->GetBinContent(iBinX, iBinY), _bgEffHist->GetBinError(iBinX, iBinY), _sigEffHist->GetBinContent(iBinX, iBinY), _sigEffHist->GetBinError(iBinX, iBinY));
		}
	}

	sbROC->SetTitle(";#epsilon_{B};#epsilon_{S}");

	return sbROC;
};


TH2* get2DCumulativeHist(TH2* _hist, Bool_t _xForward, Bool_t _yForward) {

	TH2* 												cum2dHist 	=		(TH2*)	_hist->Clone((std::string(_hist->GetName()) + "_2Dcum_" + std::to_string(_xForward) + std::to_string(_yForward)).c_str());
	cum2dHist->Reset();

	Int_t nXbins = _hist->GetNbinsX();
	Int_t nYbins = _hist->GetNbinsY();

	for (Int_t iBinX = 1; iBinX <= nXbins; iBinX++) {
		for (Int_t iBinY = 1; iBinY <= nYbins; iBinY++) {
			//////////////////////////////////////
			///////			y
			///////			|	  |
			///////			|	D | C
			///////	 		|------------------
			///////			|	A | B
			///////			|	  |
			///////			|_________________x
			///////
			Double_t iXYcum;
			Double_t iXYcumError;
			if (_xForward && _yForward) iXYcum = _hist->IntegralAndError(1, iBinX, 1, iBinY, iXYcumError);								// A
			else if (!_xForward && _yForward) iXYcum = _hist->IntegralAndError(iBinX, nXbins, 1, iBinY, iXYcumError); 					// B
			else if (_xForward && !_yForward) iXYcum = _hist->IntegralAndError(1, iBinX, iBinY, nYbins, iXYcumError);					// D
			else iXYcum = _hist->IntegralAndError(iBinX, nXbins, iBinY, nYbins, iXYcumError);											// C

			cum2dHist->SetBinContent(iBinX, iBinY, iXYcum);
			cum2dHist->SetBinError(iBinX, iBinY, iXYcumError);
		}
	}

	return cum2dHist;
};


TH1* getMaxHist(std::vector<TH1*> hists, Bool_t onlyAbs) {

	TH1* 				h1 					=		hists[0];
	TH1* 				maxHist 			= 		(TH1*) h1->Clone((std::string(h1->GetName())+"_max").c_str());
	maxHist->Reset();

	for (Int_t iBin = 1; iBin < h1->GetNbinsX(); iBin++) {

		Double_t iBinMax 					= 		- std::numeric_limits<Double_t>::max();

		for (TH1* jHist : hists) {

			Double_t iHiBinContentjHist 	= 		jHist->GetBinContent(iBin);
			if (onlyAbs) iHiBinContentjHist 	= 		std::abs(iHiBinContentjHist);
			if (iHiBinContentjHist > iBinMax) iBinMax 	= 	iHiBinContentjHist;

		}

		maxHist->SetBinContent(iBin, iBinMax);
	}

	return maxHist;
};


TH1* getInverseHist(TH1* _hist, std::string _postFix = "_Inverse") {

	TH1* inverseHist							=			(TH1*) _hist->Clone((std::string(_hist->GetName())+_postFix).c_str());
	inverseHist->Reset();

	for (Int_t iBinX = 1; iBinX <= inverseHist->GetNbinsX(); iBinX++) {
		for (Int_t iBinY = 1; iBinY <= inverseHist->GetNbinsY(); iBinY++) {
			if (_hist->GetBinContent(iBinX, iBinY) != 0) inverseHist->SetBinContent(iBinX,iBinY, 1./_hist->GetBinContent(iBinX, iBinY));
		}
	}

	return inverseHist;
};


void	drawMaxContour(TH2* _hist, Double_t _contour, std::string _contCol, Float_t _textSize, Float_t _contLineWidth, Float_t _markerSize) {
	double contoursvals[1] = {_contour};
	TH2F* 													hContour 					=		(TH2F*) _hist->Clone((std::string("cont_") + _hist->GetName() + findAndReplaceAll(removeTrailingZeros(_contour), ".", "p")).c_str());
	hContour->SetContour(1,contoursvals);
	hContour->SetFillColorAlpha(kBlack, 0.);
	hContour->Draw("CONT3 SAME");
	hContour->SetLineColor(hex2rootColor(_contCol));
	hContour->SetLineStyle(2);
	hContour->SetLineWidth(_contLineWidth);

	Int_t maxX, maxY, maxZ;
	hContour->GetMaximumBin(maxX, maxY, maxZ);

	Double_t xMax = hContour->GetXaxis()->GetBinCenter(maxX);
	Double_t yMax = hContour->GetYaxis()->GetBinCenter(maxY);

	Double_t xLeg = xMax;
	Double_t yLeg = yMax;

	Bool_t leftAlign = (xMax - hContour->GetXaxis()->GetXmin()) > (hContour->GetXaxis()->GetXmax() - xMax);
	Bool_t bottomAlign = (yMax - hContour->GetYaxis()->GetXmin()) > (hContour->GetYaxis()->GetXmax() - yMax);
	Int_t 	textAlignment = 0;
	if (leftAlign) {
		textAlignment += 30;
		xLeg -= 0.03*(hContour->GetXaxis()->GetXmax()  - hContour->GetXaxis()->GetXmin());
	} else {
		textAlignment += 10;
		xLeg += 0.03*(hContour->GetXaxis()->GetXmax()  - hContour->GetXaxis()->GetXmin());
	}
	if (bottomAlign) {
		textAlignment += 1;
		yLeg -= 0.03*(hContour->GetYaxis()->GetXmax()  - hContour->GetYaxis()->GetXmin());
	} else {
		textAlignment += 3;
		yLeg += 0.03*(hContour->GetYaxis()->GetXmax()  - hContour->GetYaxis()->GetXmin());
	}

	// TMarker	maxPoint(xMax, yMax, 22);
	// maxPoint.SetNDC();
	// maxPoint.SetMarkerSize(_markerSize);
	// maxPoint.SetMarkerColor(hex2rootColor(_contCol));
	// maxPoint.Draw("SAME");
	// maxPoint.DrawMarker(xMax, yMax);

	TLatex l;
	l.SetTextSize(_textSize);
	l.SetTextAlign(textAlignment);
	l.SetTextColor(hex2rootColor(_contCol));
	std::string 										mxStr 							=		"#bf{Max: " + ToSciString(hContour->GetBinContent(maxX, maxY), 2, "^") +
	" @ ("+to_string_with_precision(xMax,2) + ", " +to_string_with_precision(yMax,2) + ")}";
	l.DrawLatex(xLeg,yLeg,mxStr.c_str());
};


void	drawMinContour(TH2* _hist, Double_t _contour, std::string _contCol, Float_t _textSize, Float_t _contLineWidth, Float_t _markerSize) {
	double contoursvals[1] = {_contour};
	TH2F* 													hContour 					=		(TH2F*) _hist->Clone((std::string("cont_") + _hist->GetName() + findAndReplaceAll(removeTrailingZeros(_contour), ".", "p")).c_str());
	hContour->SetContour(1,contoursvals);
	hContour->SetFillColorAlpha(kBlack, 0.);
	hContour->Draw("CONT3 SAME");
	hContour->SetLineColor(hex2rootColor(_contCol));
	hContour->SetLineStyle(2);
	hContour->SetLineWidth(_contLineWidth);

	Int_t maxX, maxY, maxZ;
	hContour->GetMaximumBin(maxX, maxY, maxZ);

	Double_t xMax = hContour->GetXaxis()->GetBinCenter(maxX);
	Double_t yMax = hContour->GetYaxis()->GetBinCenter(maxY);

	Double_t xLeg = xMax;
	Double_t yLeg = yMax;

	Bool_t leftAlign = (xMax - hContour->GetXaxis()->GetXmin()) > (hContour->GetXaxis()->GetXmax() - xMax);
	Bool_t bottomAlign = (yMax - hContour->GetYaxis()->GetXmin()) > (hContour->GetYaxis()->GetXmax() - yMax);
	Int_t 	textAlignment = 0;
	if (leftAlign) {
		textAlignment += 30;
		xLeg -= 0.03*(hContour->GetXaxis()->GetXmax()  - hContour->GetXaxis()->GetXmin());
	} else {
		textAlignment += 10;
		xLeg += 0.03*(hContour->GetXaxis()->GetXmax()  - hContour->GetXaxis()->GetXmin());
	}
	if (bottomAlign) {
		textAlignment += 1;
		yLeg -= 0.03*(hContour->GetYaxis()->GetXmax()  - hContour->GetYaxis()->GetXmin());
	} else {
		textAlignment += 3;
		yLeg += 0.03*(hContour->GetYaxis()->GetXmax()  - hContour->GetYaxis()->GetXmin());
	}

	// gROOT->GetSelectedPad()->cd();
	// TMarker	maxPoint(xMax, yMax, 22);
	// maxPoint.SetMarkerSize(_markerSize);
	// maxPoint.SetMarkerColor(hex2rootColor(_contCol));
	// maxPoint.Draw("SAME");

	TLatex l;
	l.SetTextSize(_textSize);
	l.SetTextAlign(textAlignment);
	l.SetTextColor(hex2rootColor(_contCol));
	std::string 										mxStr 							=		"#bf{Max: " + ToSciString(hContour->GetBinContent(maxX, maxY), 2, "^") +
	" @ ("+to_string_with_precision(xMax,2) + ", " +to_string_with_precision(yMax,2) + ")}";
	l.DrawLatex(xLeg,yLeg,mxStr.c_str());
};


Double_t separationTMVA(const TH1* _histSig, const TH1* _histBG, Bool_t _normalize, Bool_t _preScaledByWidth) {
	////https://root.cern.ch/download/doc/tmva/TMVAUsersGuide.pdf#page=30

	if (_histSig->GetNbinsX() != _histBG->GetNbinsX()) {
		std::cout<<"Error calculating separation! Unequal bins in signal and background histograms!"<<std::endl;
		return -999;
	}

	TH1* 			sigNorm 			= 		(TH1*) _histSig->Clone((std::string(_histSig->GetName())+"_Norm").c_str());
	TH1* 			bgNorm 				= 		(TH1*) _histBG->Clone((std::string(_histBG->GetName())+"_Norm").c_str());

	if (!_preScaledByWidth) {
		sigNorm->Scale(1., "width");
		bgNorm->Scale(1., "width");
	}

	if (_normalize) {
		if (_preScaledByWidth) {

			sigNorm->Scale(1./sigNorm->Integral("width"));
			bgNorm->Scale(1./bgNorm->Integral("width"));

		} else {

			sigNorm->Scale(1./sigNorm->Integral());
			bgNorm->Scale(1./bgNorm->Integral());

		}
	}

	// std::cout<<sigNorm->Integral("width")<<"\t"<<bgNorm->Integral("width")<<std::endl;

	Double_t sigBGsep = 0.;

	for (Int_t i = 1; i <= sigNorm->GetNbinsX(); i++) {
		Double_t iNum 		= 	std::pow((sigNorm->GetBinContent(i) - bgNorm->GetBinContent(i)), 2);
		Double_t iDenom		= 	(sigNorm->GetBinContent(i) + bgNorm->GetBinContent(i));
		if (iDenom > 0.) sigBGsep += sigNorm->GetBinWidth(i) * iNum/iDenom;
	}

	sigBGsep *= 0.5;

	delete sigNorm;
	delete bgNorm;

	return sigBGsep;
};



Float_t getMaxInRange(TH1* _hist, Float_t _xMin, Float_t _xMax) {

	Float_t _maxBinVal = - std::numeric_limits<Float_t>::max();

	Float_t iStart = _hist->FindBin(_xMin);
	Float_t iEnd  = _hist->FindBin(_xMax) + 1;

	for (Int_t i = iStart; i < iEnd; i++) {
		Float_t iBinContent = _hist->GetBinContent(i);
		if (_maxBinVal < iBinContent) _maxBinVal = iBinContent;
	}

	return _maxBinVal;
};


Float_t getMinInRange(TH1* _hist, Float_t _xMin, Float_t _xMax, Float_t _zeroThres = - std::numeric_limits<Float_t>::max()) {

	Float_t _minBinVal = std::numeric_limits<Float_t>::max();

	Float_t iStart = _hist->FindBin(_xMin);
	Float_t iEnd  = _hist->FindBin(_xMax) + 1;

	for (Int_t i = iStart; i < iEnd; i++) {
		Float_t iBinContent = _hist->GetBinContent(i);
		if (iBinContent < _zeroThres) continue;
		if (_minBinVal > iBinContent) _minBinVal = iBinContent;
	}

	return _minBinVal;
};



std::pair<Double_t, Double_t>  				getHistQuantile(TH1* _hist, Double_t _quantile) {

	Double_t 			probSum[1] = {_quantile};

	Double_t 			quantileVal[1];

	_hist->GetQuantiles(1, quantileVal, probSum);

	Double_t n = _hist->GetEffectiveEntries();

	Double_t f = TMath::Gaus(quantileVal[0], _hist->GetMean(), _hist->GetStdDev(), kTRUE);

	Double_t error = 0;
	if (f > 0 && n > 1)	error = TMath::Sqrt( _quantile*(1.-_quantile)/ (n * f * f) );

	return std::make_pair(quantileVal[0], error);
};


std::pair<Double_t, Double_t>  				getHistQuantile2D(TH2* _hist, Double_t _quantile, Int_t _bin, Int_t _axis) {

	TH1D* 				projHist = nullptr;

	if (_axis == 0) {
		projHist = _hist->ProjectionX((std::string(_hist->GetName()) + "_px_" + removeTrailingZeros(rand() % 10000)).c_str(), _bin, _bin, "e");
	} else if (_axis == 1) {
		projHist = _hist->ProjectionY((std::string(_hist->GetName()) + "_pY_" + removeTrailingZeros(rand() % 10000)).c_str(), _bin, _bin, "e");
	} else {
		std::cout<<"Error! Axis must be 0 (x) or 1 (y)!"<<std::endl;

		return std::make_pair(-999, -999);
	}

	std::pair<Double_t, Double_t>  				 quantileVal = getHistQuantile(projHist, _quantile);

	delete projHist;

	return quantileVal;
};



bool TestFitSuccess(bool verbose) {
	std::string minuitstatus = std::string(gMinuit->fCstatu);
	if (minuitstatus.find("SUCCESSFUL") == std::string::npos && minuitstatus.find("STATUS=CONVERGED") == std::string::npos && minuitstatus.find("STATUS=OK") == std::string::npos) {
		if (verbose) std::cout << "  Minimization did not converge! (status_\"" << minuitstatus << "\")" << std::endl;
		return false;
	} else return true;
}



std::string randStr(UInt_t _len) {

	std::string tmp_s;
	static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	Int_t alphanumSize = sizeof(alphanum) - 1;

	srand(rand() % 1000000);

	tmp_s.reserve(_len);

	for (UInt_t i = 0; i < _len; ++i) {
		tmp_s += alphanum[rand() % (alphanumSize)];
	}

	return tmp_s;
};



std::string randHexStr(UInt_t _len) {
	std::string tmp_s;
	static const char alphanum[] = "0123456789ABCDEF";
	Int_t alphanumSize = sizeof(alphanum) - 1;

	srand(rand() % 1000000);

	tmp_s.reserve(_len);

	for (UInt_t i = 0; i < _len; ++i) {
		tmp_s += alphanum[rand() % (alphanumSize)];
	}

	return tmp_s;
};


TGraph* subtractTGraphs(TGraph* g1, TGraph* g2) {

	assert(g1->GetN() == g2->GetN());

	UInt_t nPoints = g1->GetN();

	TGraph* differenceGraph = new TGraph;

	for (UInt_t i = 0; i < nPoints; i++) {

		Double_t x1, x2, y1, y2;

		g1->GetPoint(i, x1, y1);
		g2->GetPoint(i, x2, y2);

		Double_t iDiff = 	y1 - y2;

		addPointToGraph(differenceGraph, x1, iDiff);
	}

	return differenceGraph;
};


TGraphErrors* subtractTGraphErrors(TGraphErrors* g1, TGraphErrors* g2, Bool_t g1err, Bool_t g2err) {

	assert(g1->GetN() == g2->GetN());

	UInt_t nPoints = g1->GetN();

	TGraphErrors* differenceGraph = new TGraphErrors;

	for (UInt_t i = 0; i < nPoints; i++) {

		Double_t x1, x2, y1, y2;

		g1->GetPoint(i, x1, y1);
		g2->GetPoint(i, x2, y2);

		Double_t iDiff = 	y1 - y2;

		//// square of relative error on ratio
		Double_t iDiffErr2 = 0.;

		if (g1err) iDiffErr2 += std::pow(g1->GetErrorY(i), 2);

		if (g2err) 	iDiffErr2 += std::pow(g2->GetErrorY(i), 2);

		Double_t iDiffErr = std::sqrt(iDiffErr2);

		addPointToGraph(differenceGraph, x1, g1->GetErrorX(i), iDiff, iDiffErr);
	}

	return differenceGraph;

};

TGraph* divideTGraphs(TGraph* g1, TGraph* g2) {

	assert(g1->GetN() == g2->GetN());

	UInt_t nPoints = g1->GetN();

	TGraph* ratioGraph = new TGraph;

	for (UInt_t i = 0; i < nPoints; i++) {

		Double_t x1, x2, y1, y2;

		g1->GetPoint(i, x1, y1);
		g2->GetPoint(i, x2, y2);

		Double_t iRatio = 	y1/y2;

		addPointToGraph(ratioGraph, x1, iRatio);
	}

	return ratioGraph;
};

TGraphErrors* divideTGraphErrors(TGraphErrors* g1, TGraphErrors* g2, Bool_t g1err, Bool_t g2err) {

	assert(g1->GetN() == g2->GetN());

	UInt_t nPoints = g1->GetN();

	TGraphErrors* ratioGraph = new TGraphErrors;

	for (UInt_t i = 0; i < nPoints; i++) {

		Double_t x1, x2, y1, y2;

		g1->GetPoint(i, x1, y1);
		g2->GetPoint(i, x2, y2);

		Double_t iRatio = 	y1/y2;

		//// square of relative error on ratio
		Double_t iRatioRelErr2 = 0.;

		if (g1err) iRatioRelErr2 += std::pow(g1->GetErrorY(i)/y1, 2);

		if (g2err) 	iRatioRelErr2 += std::pow(g2->GetErrorY(i)/y2, 2);

		Double_t iRatioErr = iRatio * std::sqrt(iRatioRelErr2);

		addPointToGraph(ratioGraph, x1, g1->GetErrorX(i), iRatio, iRatioErr);
	}

	return ratioGraph;
};


void divideWithoutError(TH1* _mainHist, TH1* _divisor) {

	UInt_t  													nBinsX  				=		_mainHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iMainContent  			=		_mainHist->GetBinContent(iBin);
		Double_t  												iMainError  				=		_mainHist->GetBinError(iBin);

		Double_t  												iDivisorContent  			=		_divisor->GetBinContent(iBin);

		if (iDivisorContent <= 0) {
			_mainHist->SetBinContent(iBin, -999);
			_mainHist->SetBinError(iBin, 999);
		} else {
			_mainHist->SetBinContent(iBin, iMainContent/iDivisorContent);
			_mainHist->SetBinError(iBin, iMainError/iDivisorContent);
		}
	}

};

void divideSelf(TH1* _mainHist) {

	UInt_t  													nBinsX  				=		_mainHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinContent  			=		_mainHist->GetBinContent(iBin);
		Double_t  												iBinError  				=		_mainHist->GetBinError(iBin);

		if (std::abs(iBinContent) <= 0) {

			_mainHist->SetBinContent(iBin, 0.);
			_mainHist->SetBinError(iBin, 0.);

			continue;
		}

		_mainHist->SetBinContent(iBin, 1.);
		_mainHist->SetBinError(iBin, iBinError/std::abs(iBinContent));
	}

};



void 															addContent(TH1D* _mainHist, TH1D* _addErrorsHist) {

	UInt_t  													nBinsX  				=		_mainHist->GetNbinsX();

	std::vector<Double_t>   									newContents;
	std::vector<Double_t>   									newErrors;

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinContent  			=		_mainHist->GetBinContent(iBin);
		iBinContent 																	+= _addErrorsHist->GetBinContent(iBin);

		Double_t  												iBinError  				=		_mainHist->GetBinError(iBin);

		newContents.push_back(iBinContent);
		newErrors.push_back(iBinError);
	}

	_mainHist->Reset();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		_mainHist->SetBinContent(iBin, newContents[iBin-1]);
		_mainHist->SetBinError(iBin, newErrors[iBin-1]);
	}
};


void 															addStatErrorToContent(TH1D* _variationHist, TH1D* _statErrorsHist, Float_t _addSign, Bool_t _truncateNegative) {


	UInt_t  													nBinsX  				=		_variationHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinContent  			=		_variationHist->GetBinContent(iBin) + (_addSign > 0 ? _statErrorsHist->GetBinErrorUp(iBin) : -_statErrorsHist->GetBinErrorLow(iBin));

		if (_truncateNegative) {
			(iBinContent >0 ) ? _variationHist->SetBinContent(iBin, iBinContent) : _variationHist->SetBinContent(iBin, 0.);
		} else {
			_variationHist->SetBinContent(iBin, iBinContent);
		}

	}

};


TH1D* 															addError(TH1D* _mainHist, TH1D* _addErrorsHist, std::string _newName) {

	TH1D*  														addedHist  			=		(TH1D*) _mainHist->Clone(_newName.c_str());
	addedHist->Reset();
	UInt_t  													nBinsX  				=		_mainHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinContent  			=		_mainHist->GetBinContent(iBin);
		Double_t  												iBinError  				=		std::pow(_mainHist->GetBinError(iBin),2);
		iBinError  																		+= std::pow(_addErrorsHist->GetBinError(iBin),2);
		iBinError 																		=		std::sqrt(iBinError);

		// Double_t  												iBinError  				=		_mainHist->GetBinError(iBin) + _addErrorsHist->GetBinError(iBin);

		// std::cout<<iBin<<"\t"<<iBinContent<<"\t"<<_mainHist->GetBinError(iBin)<<"\t"<<_addErrorsHist->GetBinContent(iBin)<<"\t"<<_addErrorsHist->GetBinError(iBin)<<"\t"<<iBinError<<std::endl;
		addedHist->SetBinContent(iBin, iBinContent);
		addedHist->SetBinError(iBin, iBinError);
	}

	return addedHist;
};


std::vector<TH1D*>												getPDFHessianSystSum(std::vector<TH1D*> _systHists, TH1D* _nominalHist, std::string _newName, Bool_t _addAlphaS = 1, Float_t _rescaleAlphaSunc = 1.) {

	UInt_t nHistsToSum = _systHists.size();
	if (_addAlphaS) nHistsToSum -= 2;

	TH1D*  														upHist  			=		(TH1D*) _nominalHist->Clone((_newName + "_pdfUp").c_str());
	TH1D*  														dnHist  			=		(TH1D*) _nominalHist->Clone((_newName+ "_pdfDown").c_str());
	upHist->Reset();
	dnHist->Reset();

	upHist->SetTitle((_newName+"Up").c_str());
	dnHist->SetTitle((_newName+"Dn").c_str());

	UInt_t  													nBinsX  				=		_nominalHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinNominal  			=		_nominalHist->GetBinContent(iBin);
		Double_t  												iBinError  				=		0.;

		if (iBinNominal < 0.) continue;

		for (UInt_t iSyst = 0; iSyst < nHistsToSum; iSyst++) {

			TH1D* iSystHist = _systHists[iSyst];

			Double_t  											iDeviation  			=	iSystHist->GetBinContent(iBin) - iBinNominal;
			iBinError += iDeviation*iDeviation;


		}

		iBinError = std::sqrt(iBinError);

		upHist->SetBinContent(iBin, iBinNominal + iBinError);
		dnHist->SetBinContent(iBin, iBinNominal - iBinError);

		upHist->SetBinError(iBin, _nominalHist->GetBinError(iBin));
		dnHist->SetBinError(iBin, _nominalHist->GetBinError(iBin));
	}

	if (_addAlphaS) {

		TH1D* histAlphaSdn = _systHists[nHistsToSum];
		TH1D* histAlphaSup = _systHists[nHistsToSum+1];

		TH1D* histPDFAlphaSup = (TH1D*) _systHists[nHistsToSum+1]->Clone((_newName + "_pdfAlphaSUp").c_str());;
		TH1D* histPDFAlphaSdn = (TH1D*) _systHists[nHistsToSum]->Clone((_newName + "_pdfAlphaSDown").c_str());;
		histPDFAlphaSup->Reset();
		histPDFAlphaSdn->Reset();

		for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

			Double_t  												iBinNominal  			=		_nominalHist->GetBinContent(iBin);
			Double_t  												iBinAlphaSerror			=		_rescaleAlphaSunc * std::abs(histAlphaSup->GetBinContent(iBin)  - histAlphaSdn->GetBinContent(iBin))/2.;
			Double_t  												iBinPdfError			=		std::abs(upHist->GetBinContent(iBin)  - iBinNominal);
			Double_t  												iBinPdfAlphaSerror		=		std::sqrt(iBinAlphaSerror*iBinAlphaSerror + iBinPdfError*iBinPdfError);

			histPDFAlphaSup->SetBinContent(iBin, iBinNominal + iBinPdfAlphaSerror);
			histPDFAlphaSdn->SetBinContent(iBin, iBinNominal - iBinPdfAlphaSerror);

			histPDFAlphaSup->SetBinError(iBin, _nominalHist->GetBinError(iBin));
			histPDFAlphaSdn->SetBinError(iBin, _nominalHist->GetBinError(iBin));
		}


		return std::vector<TH1D*>({upHist, dnHist, histPDFAlphaSup, histPDFAlphaSdn});
	} else {
		return std::vector<TH1D*>({upHist, dnHist});
	}
};



TH1D* 															getPDFHessianSyst(std::vector<TH1D*> _systHists, TH1D* _nominalHist, std::string _newName) {

	TH1D*  														envelopeHist  			=		(TH1D*) _nominalHist->Clone(_newName.c_str());
	envelopeHist->Reset();
	UInt_t  													nBinsX  				=		_nominalHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinNominal  			=		_nominalHist->GetBinContent(iBin);
		Double_t  												iBinError  				=		0.;

		if (iBinNominal < 0.) continue;

		for (TH1D* iSystHist : _systHists) {

			Double_t  											iDeviation  			=	iSystHist->GetBinContent(iBin) - iBinNominal;
			iBinError += iDeviation*iDeviation;


		}

		iBinError = std::sqrt(iBinError);

		// std::cout<<iBin<<"\t"<<iBinNominal<<"\t"<<iBinError<<std::endl;

		envelopeHist->SetBinContent(iBin, iBinNominal);
		envelopeHist->SetBinError(iBin, iBinError);
	}

	return envelopeHist;
};



std::vector<TH1D*>												getUpDnEnvelope(std::vector<TH1D*> _systHists, std::string _newName) {

	TH1D*  														envelopeUpHist  			=		(TH1D*) _systHists[0]->Clone((_newName+"Up").c_str());
	TH1D*  														envelopeDnHist  			=		(TH1D*) _systHists[0]->Clone((_newName+"Down").c_str());
	envelopeUpHist->Reset();
	envelopeDnHist->Reset();

	envelopeUpHist->SetTitle((_newName+"Up").c_str());
	envelopeDnHist->SetTitle((_newName+"Dn").c_str());

	UInt_t  													nBinsX  				=		envelopeUpHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinUp  				=		-std::numeric_limits<Double_t>::max();
		Double_t  												iBinDn  				=		std::numeric_limits<Double_t>::max();


		for (TH1D* iSystHist : _systHists) {

			Double_t  											iContent  				=	iSystHist->GetBinContent(iBin);
			if (iContent > iBinUp) iBinUp = iContent;
			if (iContent < iBinDn) iBinDn = iContent;
		}

		envelopeUpHist->SetBinContent(iBin, iBinUp);
		envelopeDnHist->SetBinContent(iBin, iBinDn);

		Double_t iBinErr = std::abs(iBinUp - iBinDn)/2.;

		envelopeUpHist->SetBinError(iBin, iBinErr);
		envelopeDnHist->SetBinError(iBin, iBinErr);
	}

	return std::vector<TH1D*>({envelopeUpHist, envelopeDnHist});
};


TH1D* 															getEnvelope(std::vector<TH1D*> _systHists, TH1D* _nominalHist, std::string _newName) {

	TH1D*  														envelopeHist  			=		(TH1D*) _nominalHist->Clone(_newName.c_str());
	envelopeHist->Reset();
	UInt_t  													nBinsX  				=		_nominalHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinNominal  			=		_nominalHist->GetBinContent(iBin);
		Double_t  												iBinError  				=		0.;

		if (iBinNominal <= 0.) continue;

		for (TH1D* iSystHist : _systHists) {

			Double_t  											iDeviation  			=	std::abs(iSystHist->GetBinContent(iBin) - iBinNominal);
			if (iDeviation > iBinError) iBinError = iDeviation;
		}

		envelopeHist->SetBinContent(iBin, iBinNominal);
		envelopeHist->SetBinError(iBin, iBinError);
	}

	return envelopeHist;
};



TH1D* 															getSymmetricEnvelopeError(std::vector<TH1D*> _systHists, std::string _newName) {

	TH1D*  														envelopeHist  			=		(TH1D*) _systHists[0]->Clone(_newName.c_str());
	envelopeHist->Reset();
	envelopeHist->SetTitle("");
	envelopeHist->SetName(_newName.c_str());

	UInt_t  													nBinsX  				=		_systHists[0]->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinMin = std::numeric_limits<Double_t>::max();
		Double_t  												iBinMax = -std::numeric_limits<Double_t>::max();

		for (TH1D* iSystHist : _systHists) {

			Double_t  											iSystBinConentUp  			=	iSystHist->GetBinContent(iBin) + iSystHist->GetBinError(iBin);
			Double_t  											iSystBinConentDn  			=	iSystHist->GetBinContent(iBin) - iSystHist->GetBinError(iBin);
			if (iSystBinConentUp > iBinMax) iBinMax = iSystBinConentUp;
			if (iSystBinConentDn < iBinMin) iBinMin = iSystBinConentDn;
		}

		Double_t  												iBinContent = (iBinMax + iBinMin)/2.;
		Double_t  												iBinErr = (iBinMax - iBinMin)/2.;

		if (iBinContent <= 0) continue;

		envelopeHist->SetBinContent(iBin, iBinContent);
		envelopeHist->SetBinError(iBin, iBinErr);
	}

	return envelopeHist;
};


TH1D* 															getSymmetricEnvelope(std::vector<TH1D*> _systHists, std::string _newName) {

	TH1D*  														envelopeHist  			=		(TH1D*) _systHists[0]->Clone(_newName.c_str());
	envelopeHist->Reset();
	envelopeHist->SetTitle("");
	envelopeHist->SetName(_newName.c_str());

	UInt_t  													nBinsX  				=		_systHists[0]->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinMin = std::numeric_limits<Double_t>::max();
		Double_t  												iBinMax = -std::numeric_limits<Double_t>::max();

		for (TH1D* iSystHist : _systHists) {

			Double_t  											iSystBinConent  			=	iSystHist->GetBinContent(iBin);
			if (iSystBinConent > iBinMax) iBinMax = iSystBinConent;
			if (iSystBinConent < iBinMin) iBinMin = iSystBinConent;
		}

		Double_t  												iBinContent = (iBinMax + iBinMin)/2.;
		Double_t  												iBinErr = (iBinMax - iBinMin)/2.;

		if (iBinContent <= 0) continue;

		envelopeHist->SetBinContent(iBin, iBinContent);
		envelopeHist->SetBinError(iBin, iBinErr);
	}

	return envelopeHist;
};


TH1D* 															addEnvelope(std::vector<TH1D*> _systHists, TH1D* _nominalHist, TH1D* _prevEnvelop, std::string _newName, Bool_t _addQuad) {

	TH1D*  														newEnvelopeHist 		=		(TH1D*) _systHists[0]->Clone(_newName.c_str());
	newEnvelopeHist->Reset();
	newEnvelopeHist->SetTitle("");
	newEnvelopeHist->SetName(_newName.c_str());

	UInt_t  													nBinsX  				=		_nominalHist->GetNbinsX();

	for (UInt_t iBin = 1; iBin <= nBinsX; iBin++) {

		Double_t  												iBinMin = std::numeric_limits<Double_t>::max();
		Double_t  												iBinMax = -std::numeric_limits<Double_t>::max();

		for (TH1D* iSystHist : _systHists) {

			Double_t  											iSystBinConent  			=	iSystHist->GetBinContent(iBin);
			if (iSystBinConent > iBinMax) iBinMax = iSystBinConent;
			if (iSystBinConent < iBinMin) iBinMin = iSystBinConent;
		}

		Double_t systErrUp = iBinMax - _nominalHist->GetBinContent(iBin);
		Double_t systErrDn = -(iBinMin - _nominalHist->GetBinContent(iBin));

		if (_addQuad) {

			Double_t combinedErrUp = _prevEnvelop->GetBinContent(iBin) + _prevEnvelop->GetBinError(iBin) - _nominalHist->GetBinContent(iBin);
			Double_t combinedErrDn = -(_prevEnvelop->GetBinContent(iBin) - _prevEnvelop->GetBinError(iBin) - _nominalHist->GetBinContent(iBin));

			if (systErrUp > 0) {
				combinedErrUp = std::sqrt(combinedErrUp*combinedErrUp + systErrUp*systErrUp);

			} else {
				combinedErrDn = std::sqrt(combinedErrDn*combinedErrDn + systErrUp*systErrUp);
			}


			if (systErrDn > 0) {
				combinedErrDn = std::sqrt(combinedErrDn*combinedErrDn + systErrDn*systErrDn);

			} else {
				combinedErrDn = std::sqrt(combinedErrUp*combinedErrUp + systErrDn*systErrDn);
			}

			iBinMin += combinedErrUp;
			iBinMin -= combinedErrDn;

		} else {

			iBinMin = _prevEnvelop->GetBinContent(iBin) - _prevEnvelop->GetBinError(iBin);
			iBinMax = _prevEnvelop->GetBinContent(iBin) + _prevEnvelop->GetBinError(iBin);

			if (systErrDn > 0) {
				iBinMin -= systErrDn;
			} else {
				iBinMax += systErrDn;
			}

			if (systErrUp > 0) {
				iBinMax += systErrUp;
			} else {
				iBinMin -= systErrUp;
			}


		}


		Double_t  												iBinErr = (iBinMax - iBinMin)/2.;
		Double_t  												iBinContent = (iBinMax + iBinMin)/2.;


		newEnvelopeHist->SetBinContent(iBin, iBinContent);
		newEnvelopeHist->SetBinError(iBin, iBinErr);
	}

	return newEnvelopeHist;
};


template <class HISTCLASS>
std::vector<Double_t> getHistMinMax(std::vector<HISTCLASS*> _hists, Bool_t _ignoreZero, Bool_t _ignoreError) {

	Double_t  												yMin = std::numeric_limits<Double_t>::max();
	Double_t  												yMax = -std::numeric_limits<Double_t>::max();

	for (Int_t iBin = 1; iBin <= _hists[0]->GetNbinsX(); iBin++) {

		for (TH1* iHist : _hists) {

			Double_t  											iBinUp  			=	iHist->GetBinContent(iBin) + (_ignoreError ? 0. : iHist->GetBinError(iBin));
			Double_t  											iBinDn  			=	iHist->GetBinContent(iBin) - (_ignoreError ? 0. : iHist->GetBinError(iBin));

			if (iBinUp > yMax && (!_ignoreZero || iBinUp > 0)) yMax = iBinUp;
			if (iBinDn < yMin && (!_ignoreZero || iBinDn > 0)) yMin = iBinDn;
		}


	}

	return std::vector<Double_t>({yMin, yMax});
};



void scaleHistErrBandByTanh(TH1* _hist, Double_t _yCentre) {

	for (Int_t i = 1; i < _hist->GetNbinsX(); i++) {

		Double_t iBinUp = _hist->GetBinContent(i) + _hist->GetBinError(i);
		Double_t iBinDn = _hist->GetBinContent(i) + _hist->GetBinError(i);

		/// new up/dn/central
		iBinUp = tanh(iBinUp - _yCentre);
		iBinDn = tanh(_yCentre - iBinDn);
		Double_t iBinContent = (iBinUp + iBinDn)/2.;
		Double_t iBinErr = (iBinUp - iBinDn)/2.;

		_hist->SetBinContent(i, iBinContent);
		_hist->SetBinError(i, iBinErr);
	}

};


TGaxis* makeTanhYaxis(TPad* _pad, Double_t _yMin, Double_t _yMax, Double_t _yCentre) {

	std::string fName = randStr();
	std::string fString = "tanh(x - " + std::to_string(_yCentre) + ")+" + std::to_string(_yCentre);
	TF1* tanhAxisFunc = new TF1(fName.c_str(), fString.c_str(), _yMin, _yMax);

	TGaxis *axis = new TGaxis(_pad->GetUxmin(), _pad->GetUymin(), _pad->GetUxmin(), _pad->GetUymax(), fName.c_str(), 505);

	return axis;
};




void moveOverflowToLastBin(TH1* _hist) {

	Double_t overflowBinContent  = _hist->GetBinContent(_hist->GetNbinsX()+1);
	Double_t overflowBinError    = _hist->GetBinContent(_hist->GetNbinsX()+1);

	Double_t lastBinContent      = _hist->GetBinContent(_hist->GetNbinsX());
	Double_t lastBinError        = _hist->GetBinContent(_hist->GetNbinsX());

	Double_t newLastBinContent 	 = overflowBinContent + lastBinContent;
	Double_t newLastBinError 	   = std::sqrt(overflowBinError*overflowBinError + lastBinError+lastBinError);

	_hist->SetBinContent(_hist->GetNbinsX()+1, 0.);
	_hist->SetBinError(_hist->GetNbinsX()+1, 0.);

	_hist->SetBinContent(_hist->GetNbinsX(), newLastBinContent);
	_hist->SetBinError(_hist->GetNbinsX(), newLastBinError);
};


void moveUnderflowToFirstBin(TH1* _hist) {

	Double_t underflowBinContent  = _hist->GetBinContent(0);
	Double_t underflowBinError    = _hist->GetBinContent(0);

	Double_t firstBinContent      = _hist->GetBinContent(1);
	Double_t firstBinError        = _hist->GetBinContent(1);

	Double_t newFirstBinContent 	 = underflowBinContent + firstBinContent;
	Double_t newFirstBinError 	   = std::sqrt(underflowBinError*underflowBinError + firstBinError+firstBinError);

	_hist->SetBinContent(0, 0.);
	_hist->SetBinError(0, 0.);

	_hist->SetBinContent(1, newFirstBinContent);
	_hist->SetBinError(1, newFirstBinError);

};


Double_t getWeightedBinVariation(std::vector<TH1*> _variations, TH1* _nominalHist) {

	Double_t weightedVariation = 0.;

	Double_t iNomIntegral = _nominalHist->Integral();

	for (Int_t i = 1; i <= _nominalHist->GetNbinsX(); i++) {

		Double_t iNomBinContent = _nominalHist->GetBinContent(i);
		Double_t iVariation 	= 0.;


		for (TH1* jVarHist : _variations) {

			iVariation += std::abs(jVarHist->Integral() - iNomIntegral);

			// iVariation += std::abs(jVarHist->GetBinContent(i) - iNomBinContent)/iNomBinContent;
		}

		weightedVariation += iVariation;
	}

	return weightedVariation/iNomIntegral;
};


Double_t getHistSupremum(TH1* _hist, const Double_t _upperBound) {

	Double_t hMax = -std::numeric_limits<Double_t>::max();

	for (Int_t i = 1; i <= _hist->GetNbinsX(); i++) {
		Double_t iContent = _hist->GetBinContent(i);

		if (iContent > _upperBound) continue;

		if (iContent > hMax) hMax = iContent;
	}

	return hMax;
};


Double_t getHistInfimum(TH1* _hist, const Double_t _lowerBound) {

	Double_t hMin = std::numeric_limits<Double_t>::max();

	for (Int_t i = 1; i <= _hist->GetNbinsX(); i++) {
		Double_t iContent = _hist->GetBinContent(i);

		if (iContent < _lowerBound) continue;

		if (iContent < hMin) hMin = iContent;
	}

	return hMin;
};


void deleteObjectFromTFile(std::string _filePath, std::string _objName) {

	TFile f(_filePath.c_str(), "UPDATE");
	f.Delete((_objName + ";*").c_str());
	f.Write();
	f.Close();

};


void copyFile(std::string _sourcePath, std::string _destPath) {

	std::ifstream src(_sourcePath, std::ios::binary);
	std::ofstream dest(_destPath, std::ios::binary);
	dest << src.rdbuf();

};


void writeStringToFile(std::string _stringToWrite, std::string _filePath, Bool_t _append, Bool_t _verbose) {
	std::ofstream outFile;
	if (_append && file_exists(_filePath)) {
		outFile.open(_filePath, std::ios_base::app);
	} else {
		outFile.open(_filePath, std::ios_base::out);
	}
	outFile << _stringToWrite;
	outFile.close();

	if (_verbose) std::cout<<"Wrote string to file ("<<_filePath<<") in append mode = "<<_append<<std::endl;
};


TGraphAsymmErrors* getPoissonErrorRatioGraph(TH1* _dataHist, TH1* _predHist, Float_t _confidence){
	Double_t                                                  alphaOver2 = (1. - _confidence)/2.;
	TGraphAsymmErrors* dataRatioPoisErrPlot = new TGraphAsymmErrors();
	for (Int_t i = 1 ; i <= _dataHist->GetNbinsX(); i++) {
		if (_predHist->GetBinContent(i) <= 0) continue;
		Double_t iYUp;
		Double_t iYDn;
		Double_t iY = _dataHist->GetBinContent(i);
		if (iY > 0) {
                ////  Casella-Berger's Statistical Inference : Example 9.2.15,
			iYUp      = ROOT::Math::gamma_quantile_c(alphaOver2, iY+1., 1.);
			iYDn      = ROOT::Math::gamma_quantile(alphaOver2, iY, 1.);
		} else continue;
		iY        /= _predHist->GetBinContent(i);
		iYUp      /= _predHist->GetBinContent(i);
		iYDn      /= _predHist->GetBinContent(i);
		Double_t iX        = _dataHist->GetBinCenter(i);
		Double_t iXErr     = _dataHist->GetBinWidth(i)/2.;
		Double_t iYErrUp   = iYUp - iY;
		Double_t iYErrDn   = iY - iYDn;
		addPointToGraph(dataRatioPoisErrPlot, iX, iXErr, iXErr, iY, iYErrDn, iYErrUp);
	}
	return dataRatioPoisErrPlot;
}


TGraphAsymmErrors* getPoissonErrorGraph(TH1* _hist, Bool_t _divideByWidth, Float_t _confidence) {
	////https://people.maths.ox.ac.uk/gilesm/talks/poisson_2013.pdf
	////https://math.stackexchange.com/questions/467341/question-about-connection-between-poisson-and-gamma-distributions
	////Casella-Berger's Statistical Inference : Example 9.2.15,
	Double_t alpha = (1. - _confidence)/2.;
	TGraphAsymmErrors* graph = new TGraphAsymmErrors();
	for (Int_t i = 1; i <= _hist->GetNbinsX(); i++) {
		Double_t iYUp;
		Double_t iYDn;
		Double_t iY = _hist->GetBinContent(i);
		if (iY >= 0) {
			iYUp      = ROOT::Math::gamma_quantile_c(alpha, iY+1., 1.);
			iYDn      = ROOT::Math::gamma_quantile(alpha, iY, 1.);
		} else continue;

		Double_t iBinWidth = _hist->GetBinWidth(i);
		Double_t iX 	   = _hist->GetBinCenter(i);
		if (_divideByWidth) {
			iY /= iBinWidth;
			iYUp /= iBinWidth;
			iYDn /= iBinWidth;
		}
		Double_t iYErrUp   = iYUp - iY;
		Double_t iYErrDn   = iY - iYDn;
		Double_t iXErr     = iBinWidth/2.;
		addPointToGraph(graph, iX, iXErr, iXErr, iY, iYErrDn, iYErrUp);
	}
	return graph;
};


template <typename T>
Bool_t areFloatsEqual(T val1, T val2, T tmpEpsilon) {

	Bool_t tmpFloatsAreEqual = (std::abs(val1 - val2) < tmpEpsilon);
	return tmpFloatsAreEqual;
};


template <typename T>
Bool_t areVectorsEqual(std::vector<T> vec1, std::vector<T> vec2, T tmpEpsilon) {

	Bool_t tmpVectorsAreEqual = (vec1.size() == vec2.size());
	if (!tmpVectorsAreEqual) return tmpVectorsAreEqual;

	for (UInt_t i = 0; i < vec1.size(); i++) {
		tmpVectorsAreEqual = (tmpVectorsAreEqual && areFloatsEqual(vec1[i], vec2[i], tmpEpsilon));
		if (!tmpVectorsAreEqual) break;
	}

	return tmpVectorsAreEqual;
};

#endif
