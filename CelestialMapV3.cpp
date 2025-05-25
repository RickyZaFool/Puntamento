#include <iostream>
#include <fstream>
#include <TGraphPolar.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphPolargram.h>
#include <TMath.h>
#include <TLatex.h>
#include <TStyle.h>
#include <cmath>
#include <ctime>

//Useful variables and conversions
auto pi = TMath::Pi();
double degToRad = pi / 180;


//JDN calculation
long JDN(int y, int m, int d) {
    if (m <= 2) {
        y -= 1;
        m += 12;
    }
    return (1461 * (y + 4800)) / 4
         + (367 * (m - 2)) / 12
         - (3 * ((y + 4900) / 100)) / 4
         + d - 32075;
}

bool is_last_sunday(int year, int month, int day) {
    std::tm timeinfo = {};
    timeinfo.tm_year = year - 1900;
    timeinfo.tm_mon = month - 1;
    timeinfo.tm_mday = day;

    std::mktime(&timeinfo); // Normalize timeinfo to get weekday
    int weekday = timeinfo.tm_wday;

    // Find the last day of the month
    int last_day;
    if (month == 2) {
        // Leap year check
        bool leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        last_day = leap ? 29 : 28;
    } else if (month == 4 || month == 6 || month == 9 || month == 11) {
        last_day = 30;
    } else {
        last_day = 31;
    }

    return weekday == 0 && day + 7 > last_day;
}

bool is_dst_italy(int year, int month, int day, int hour) {
    // Find DST start: last Sunday of March at 02:00
    // Find DST end: last Sunday of October at 03:00

    if (month < 3 || month > 10) return false;
    if (month > 3 && month < 10) return true;

    if (month == 3) {
        // March: check if it's after the last Sunday at 02:00
        for (int d = 31; d >= 25; --d) {
            if (is_last_sunday(year, 3, d)) {
                if (day > d) return true;
                if (day == d && hour >= 2) return true;
                return false;
            }
        }
    }

    if (month == 10) {
        // October: check if it's before the last Sunday at 03:00
        for (int d = 31; d >= 25; --d) {
            if (is_last_sunday(year, 10, d)) {
                if (day < d) return true;
                if (day == d && hour < 3) return true;
                return false;
            }
        }
    }

    return false;
}


double local_sidereal_time(double jd, double longitude) {
    // Step 1: Julian Century (T)
    double T = (jd - 2451545.0) / 36525.0;

    // Step 2: Greenwich Mean Sidereal Time (GMST) in degrees
    double GMST = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 
                  T * T * 0.000387933 - T * T * T / 38710000;
    
    // Make sure GMST is within 0 to 360 degrees
    GMST = fmod(GMST, 360.0);

    // Step 3: Convert GMST to Local Sidereal Time (LST) in degrees
    double LST = GMST + longitude;
    
    // Step 4: Adjust LST to be within 0 to 360 degrees
    LST = fmod(LST, 360.0);
    if (LST < 0) {
        LST += 360.0;
    }

    return LST;
}

double calculateAltitude(double latitude, double declination, double hourAngle){
    double alt = TMath::ASin(TMath::Sin(latitude * degToRad)*TMath::Sin(declination * degToRad) + TMath::Cos(latitude * degToRad)*TMath::Cos(declination * degToRad)*TMath::Cos(hourAngle * degToRad));
    return alt / degToRad;
} 

double calculateAzimuth(double declination, double hourAngle, double altitude, double latitude){
    double sinAz = -TMath::Sin(hourAngle * degToRad) * TMath::Cos(declination * degToRad);
    double cosAz = -(TMath::Cos(hourAngle * degToRad) * TMath::Cos(declination * degToRad) * TMath::Sin(latitude * degToRad) - TMath::Sin(declination * degToRad) * TMath::Cos(latitude * degToRad));
    double azimuth = atan2(sinAz, cosAz);

    return azimuth;
}

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
    TLatex *LabelStellariumTitle = new TLatex (0.128, 1.05, "Stellarium");
    LabelStellariumTitle->SetTextAlign(25);
    LabelStellariumTitle->SetTextSize(0.02);
    LabelStellariumTitle->Draw();
    TLatex *LabelObservatoryTitle = new TLatex ( - 0.128, 1.05, "Obserbvatory");
    LabelObservatoryTitle->SetTextAlign(25);
    LabelObservatoryTitle->SetTextSize(0.02);
    LabelObservatoryTitle->Draw();
    for(int i = 0; i<RadialDivs; i++){
        double y = (i+1) * 0.1 + 0.0125*i; //TODO: need better alignment, the division is non linear.
        double xStellarium = 0.02 + 0.012 * i;
        TLatex *LabelStellarium = new TLatex(xStellarium, y, std::to_string(90 - (i+1)*10).c_str());
        LabelStellarium->SetTextAlign(25);
        LabelStellarium->SetTextSize(0.02);
        LabelStellarium->Draw();
        double xObservatory = - (0.02 + 0.012 * i);
        TLatex *LabelObservatory = new TLatex(xObservatory, y, std::to_string((i+1)*10).c_str());
        LabelObservatory->SetTextAlign(25);
        LabelObservatory->SetTextSize(0.02);
        LabelObservatory->Draw();
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
    
    double RaH = 0, RaM = 0, RaS = 0, DecD = 0, DecM = 0, DecS = 0, startTime = 0, endTime = 0;
    int day = 0, month = 0, year = 0;
    std::vector<double> raList(0);
    std::vector<double> decList(0);
    std::vector<double> startTimeList(0);
    std::vector<double> endTimeList(0);
    predictedFile >> day >> month >> year;
    
    long jdn = JDN(year, month, day);
    std::vector<double> JDList;


    
    while(predictedFile >> RaH >> RaM >> RaS >> DecD >> DecM >> DecS >> startTime >> endTime ){
    	double raInDec = (RaH + RaM / 60 + RaS / 3600)*15;
    	double decInDec = DecD + DecM / 60 + DecS / 3600;
    	raList.push_back(raInDec);
    	decList.push_back(decInDec);
    	double startHour = std::floor(startTime);
    	double startMinutes = (startTime - startHour) * 100;
    	double endHour = std::floor(endTime);
    	double endMinutes = (endTime - endHour) * 100;
        double currJdn = jdn;
        int legal;
        if(is_dst_italy(year, month, day, startHour)){
            legal = -2;
        }else{
            legal = -1;
        }
        startHour += legal;
        if(startHour < 0){
            startHour += 24;
            currJdn -= 1;
        }
        if(is_dst_italy(year, month, day, endHour)){
            legal = -2;
        }else{
            legal = -1;
        }
        endHour += legal;
        if(endHour < 0){
            endHour += 24;
        }

        double JD = currJdn + (startHour - 12) / 24 + startMinutes / 1440;
        JDList.push_back(JD);
    	startTimeList.push_back(startHour * 3600 + startMinutes  * 60);
		endTimeList.push_back(endHour * 3600 + endMinutes * 60);
    }
    
    
    int nSteps = 1000;

    TGraphPolar grs[startTimeList.size()];
    
    for(unsigned long int i=0; i<startTimeList.size(); i++){

        std::vector<double> AzPredictedList(0);
        std::vector<double> ZdPredictedList(0);       

        if(endTimeList[i] - startTimeList[i] < 0){
            endTimeList[i] += 24 * 3600;
        }
        std::cout << endTimeList[i] << " " << startTimeList[i] << std::endl;
        std::cout << (endTimeList[i] - startTimeList[i])/3600 << std::endl;

    	double LST = local_sidereal_time(JDList[i],ObservatoryLongitude);
        double hourAngle = LST - raList[i];
        hourAngle = fmod(hourAngle, 360.0);
        double altitude = calculateAltitude(CelestialPoleX , decList[i], hourAngle);
        ZdPredictedList.push_back(90 - altitude);
        double azimuthPredicted = calculateAzimuth(decList[i], hourAngle, altitude, CelestialPoleX);
        AzPredictedList.push_back(azimuthPredicted  + directionOfNorth * degToRad);
        for(int j = 0; j < nSteps; j++){
    		double timePassed = (endTimeList[i] - startTimeList[i]) * double(j + 1) / nSteps;
            LST = local_sidereal_time(JDList[i] + timePassed / 86400.0 , ObservatoryLongitude);
            hourAngle = LST - raList[i];
            altitude = calculateAltitude(CelestialPoleX, decList[i], hourAngle);
            ZdPredictedList.push_back(90 - altitude);
            azimuthPredicted = calculateAzimuth(decList[i], hourAngle, altitude, CelestialPoleX);
            AzPredictedList.push_back(azimuthPredicted  + directionOfNorth * degToRad);
    	}
        Double_t AzPredictedGraph[AzPredictedList.size()];  //Thanks, root
        Double_t ZdPredictedGraph[ZdPredictedList.size()];  //Thanks, root

        for(unsigned long int i = 0; i < AzPredictedList.size(); i++){
            AzPredictedGraph[i] = AzPredictedList[i];
            ZdPredictedGraph[i] = ZdPredictedList[i];
        }
        grs[i] = TGraphPolar(AzPredictedList.size(), AzPredictedGraph, ZdPredictedGraph);
        grs[i].SetPolargram(polarFrame);
        grs[i].SetLineColor(kRed);
        grs[i].SetMinRadial(0);
        grs[i].SetMaxRadial(90);
        grs[i].SetLineWidth(1.7);
        grs[i].Draw("L SAME");
        c1.Update();
    }   

    //ROOT graphing tool. THESE LINES MUST BE AT THE END OF THE PAGE ALWAYS.
    
    c1.Update();
    app.Run();

    return 0;
}
