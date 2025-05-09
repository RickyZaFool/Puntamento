#include <iostream>
#include <fstream>
#include <TGraphPolar.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphPolargram.h>
#include <TMath.h>
#include <TLatex.h>
#include <TStyle.h>

std::string LabelsPolarTextual24Divs[24] = {
    "N",
    "15",
    "30",
    "NE",
    "60",
    "75",
    "E",
    "105",
    "120",
    "SE",
    "150",
    "165",
    "S",
    "195",
    "210",
    "SW",
    "240",
    "255",
    "W",
    "285",
    "300",
    "NW",
    "330",
    "345"
};

int main(){
	//Useful variables and conversions
	auto pi = TMath::Pi();
    double degToRad = pi / 180;
    
    //NOTE! The ROOT standard for polar graphing has north to the right.
    //Thus, the Azimuth 0 is parallel to positive X in ROOT's internal system.
    //Thanks, ROOT.
    double ObservatoryLongitude = 9.1986191; //(thanks, google maps) :)
    double CelestialPoleY = 0; //  by definition
    double CelestialPoleX = 44.5892595; //Latitude of the observatory (thanks, google maps) :)
    
    //Root application initialization
    TApplication app("app", 0,0);   //initialization
    gStyle->SetOptTitle(0);			//removing title
    TCanvas c1;						//the canvas for graphing
    
    //Graph style choices
    auto polarFrame = new TGraphPolargram("polarFrame", 0, 360*degToRad, 0, 90);
    polarFrame->SetToRadian();
    polarFrame->SetNdivRadial(9);
    polarFrame->SetNdivPolar(24); //STANDARD IS 24. ANY OTHER DIV NUMBER WILL NOT HAVE LETTER LABELING
    polarFrame->Draw("P SAME");
    polarFrame->SetPolarLabelSize(0);
    double r = 1.1; // Relative distance from graph center, used to label the divisions.
    double directionOfNorth = 90; //Root's default is right.
    int PolarDivs = polarFrame->GetNdivPolar();
    int RadialDivs = polarFrame->GetNdivRadial();
    double divAngle = 360.0/PolarDivs; // Angle of a single polar division
    
    //Labelling of polar divisions.
    if(PolarDivs == 24){
	    for(int i = 0; i<PolarDivs; i++){
    	    double x = r * TMath::Cos((directionOfNorth + i * divAngle) * degToRad);
    	    double y = r * TMath::Sin((directionOfNorth + i * divAngle) * degToRad);
    	    TLatex *Label = new TLatex(x, y, LabelsPolarTextual24Divs[i].c_str());
    	    Label->SetTextAlign(22); // Center align
    	    Label->SetTextSize(0.05);
    	    Label->Draw();
    	}
    }
    else{
	    for(int i = 0; i<PolarDivs; i++){
    	    double x = r * TMath::Cos((directionOfNorth + i * divAngle) * degToRad);
    	    double y = r * TMath::Sin((directionOfNorth + i * divAngle) * degToRad);
    	    TLatex *Label = new TLatex(x, y, std::to_string(i).c_str());
    	    Label->SetTextAlign(22); // Center align
    	    Label->SetTextSize(0.05);
    	    Label->Draw();
    	}
    }
    
    //Labelling of radial division. 
    //In the telescope system, 90 is the horizon. In stellarium, 0 is the orizon. TODO: Display both with corresponding labelling 
    for(int i = 0; i<RadialDivs; i++){
        double y = (i+1) * 0.1 + 0.0125*i; //TODO: need better alignment, the division is non linear.
        double x = 0;
        TLatex *Label = new TLatex(x, y, std::to_string((i+1)*10).c_str());
        Label->SetTextAlign(25);
        Label->SetTextSize(0.02);
        Label->Draw();
    }
    
    //File handling
    
    std::ifstream takenFile("taken.dat");
    std::ifstream predictedFile("predicted.dat");
    
    double takenAz = 0, takenZd = 0;
    
    std::vector<double> takenAzList(0);
    std::vector<double> takenZdList(0);
    
    //Taken stars processing start
    
    while(takenFile >> takenAz >> takenZd){
    	takenAzList.push_back((takenAz + directionOfNorth)*degToRad); //Thanks, root. You really had to default north right.
    	takenZdList.push_back(takenZd);
    }
    
    //Brief explanation of the following few lines: ROOT is strange in its way of handling graphs (albeit really useful)
    //The _t type of variables is ROOT internal stuff, and its constructors work way better when using that type of variables.
    //Whenever a conversion is made that seems strange or useless, it's done to work around constructor limitations
    
    Double_t takenAzGraph[takenAzList.size()];  //Thanks, root
    Double_t takenZdGraph[takenZdList.size()];  //Thanks, root
    
    for(unsigned long int i=0; i<takenAzList.size(); i++){
        takenAzGraph[i] = takenAzList[i];
        takenZdGraph[i] = takenZdList[i];
    }
    
    //Graph creation and styling
    TGraphPolar gr1(takenAzList.size(), takenAzGraph, takenZdGraph);
    
    gr1.SetPolargram(polarFrame);
    gr1.SetMarkerStyle(29);
    gr1.SetMarkerColor(kBlack);
    gr1.SetMinRadial(0);
    gr1.SetMaxRadial(90);
    gr1.SetMarkerSize(1.5);
    gr1.Draw("P SAME");
    
    //Taken stars processing end
    
    //Prediction processing
    
    double RaH = 0, RaM = 0, RaS = 0, DecD = 0, DecM = 0, DecS = 0, startTime = 0, endTime = 0, day = 0, month = 0, year = 0;
    vector<double> raList(0);
    vector<double> decList(0);
    vector<double> startTimeList(0);
    vector<double> endTimeList(0);
    predictedFile >> day >> month >> year;
    
    double JDN = std::floor(1461*(year + 4800 + stf::floor((month - 14)/12))/4) + std::floor(367*(month - 2 - 12 * std::floor((month-14)/12))/12) - std::floor(3*(std::floor((year + 4900 + std::floor((month-14)/12))/100))/4) + day - 32075;
    
    double legal;
    if(){
    	legal = 1;
    }
    else{
    	legal = 2;
    }
    
    double UT = hour + legal;
    
    double JD = 
    
    while(predictedFile >> Rah >> RaM >> RaS >> DecH >> DecM >> DecS >> startTime >> endTime ){
    	double raInDec = (Rah + RaM / 60 + RaS / 3600)*15;
    	double decInDec = DecD + DecM / 60 + DecS / 3600;
    	raList.push_back(raInDec);
    	decList.push_back(decInDec);
    	double startHour = std::floor(startTime);
    	double startMinutes = (startTime - startHour) * 100;
    	double endHour = std::floor(endTime);
    	double endMinutes = (endTime - endHour) * 100;
    	startTimeList.push_back(startHour * 3600 + startMinutes  * 60);
		endTimeList.push_back(endHour * 3600 + endMinutes * 60);
    }
    
    int nSteps = 10000;
    
    for(unsigned long int i=0; i<startTimeList.size(); i++){
    	for(int j = 0; j < nSteps; j++){
    		
    	}
    }
    
    
    
    
    
    
    
    
    Double_t AzPredictedGraph[];  //Thanks, root
    Double_t ZdPredictedGraph[];  //Thanks, root



    
    //Graph of the start and end points

    TGraphPolar gr2(predictionAzList.size(), AzPredictedGraph, ZdPredictedGraph);
    gr2.SetPolargram(polarFrame);
    gr2.SetMarkerStyle(29);
    gr2.SetMarkerColor(kRed);
    gr2.SetMinRadial(0);
    gr2.SetMaxRadial(90);
    gr2.SetMarkerSize(1.5);
    gr2.Draw("P SAME");
    
    
    
    
    
    
    
    //ROOT graphing tool. THESE LINES MUST BE AT THE END OF THE PAGE ALWAYS.
    
    c1.Update();
    app.Run();

    return 0;
}
