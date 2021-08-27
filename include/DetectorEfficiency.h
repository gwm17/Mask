#ifndef DETECTOREFFICIENCY_H
#define DETECTOREFFICIENCY_H

#include <THashTable.h>
#include <TH1.h>
#include <TH2.h>
#include <string>

class DetectorEfficiency {
public:
	DetectorEfficiency() { m_rxn_type = -1; };
	virtual ~DetectorEfficiency() {};

	inline void SetReactionType(int rxntype) { m_rxn_type = rxntype; };
	virtual void CalculateEfficiency(const std::string& filename) = 0;
	virtual void DrawDetectorSystem(const std::string& filename) = 0;
	virtual double RunConsistencyCheck() = 0;

protected:
	void MyFill(THashTable* table, const std::string& name, const std::string& title, int bins, float min, float max, double val) {
		TH1F* h = (TH1F*) table->FindObject(name.c_str());
		if(h) {
			h->Fill(val);
		} else {
			h = new TH1F(name.c_str(), title.c_str(), bins, min, max);
			h->Fill(val);
			table->Add(h);
		}
	}
    void MyFill(THashTable* table, const std::string& name, const std::string& title, int binsx, float minx, float maxx, double valx, int binsy, float miny, float maxy, double valy) {
    	TH2F* h = (TH2F*) table->FindObject(name.c_str());
		if(h) {
			h->Fill(valx, valy);
		} else {
			h = new TH2F(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
			h->Fill(valx, valy);
			table->Add(h);
		}
    }
	inline bool IsDoubleEqual(double x, double y) { return std::fabs(x-y) < epsilon ? true : false; };


    virtual void Run2Step(const std::string& filename) = 0;
	virtual void Run3Step(const std::string& filename) = 0;
	virtual void RunDecay(const std::string& filename) = 0;
	virtual void Run1Step(const std::string& filename) = 0;

	static constexpr double epsilon = 1.0e-6;
	int m_rxn_type;
};

#endif