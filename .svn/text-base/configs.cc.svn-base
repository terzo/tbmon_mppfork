
/*void DESY2011(TbConfig &config){

  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  eut->addMatchDUT(10);
  eut->addMatchDUT(20);
  eut->nMatches(1);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  DUT* m10 = makeFEI4("FEI4", 10, 3);
  m10->lv1Range(0,10);
  config.addDut(m10);
  DUT* m20 = makeFEI4("FEI4", 20, 3);
  m20->lv1Range(0,10);
  config.addDut(m20);
}*/
void DESY2011(TbConfig &config){
  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  Translator* translator = new Translator();
  //eut->addTranslation("calibs/run30080-30088-checktranslations-13-translation.txt",13,0,0);
  //eut->addTranslation("calibs/run30080-30088-checktranslations-14-translation.txt",14,0,0);
  eut->addMatchDUT(13);
  eut->addMatchDUT(14);
  eut->nMatches(1);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  DUT* m13 = makeFEI3("FEI3", 13, 3);
  m13->lv1Range(0,10);
  config.addDut(m13);
  DUT* m14 = makeFEI4("FEI4", 14, 3);
  m14->lv1Range(0,10);
  config.addDut(m14);
  }

void BAT2008(TbConfig &config){
#ifdef SITEBAT2008_SET
  siteBAT2008(config);
#else 
  cout << "Site config for Bat2008 not found" << endl;
#endif
  config.addBuilder(new BatTrack);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder());
  config.addBuilder(new CheckRegion);

  Chi2Builder* chi2 = new Chi2Builder(15.0); //Chi2 cut at 20.
  config.addBuilder( (EventBuilder*) chi2);
  DUT* m160 = new DUT("PLANAR", 160,0);
  //m162->lv1Range(3,10);
  //m162->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-162-masks.txt");
  //m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-162-etacalib.txt");
  config.addDut(m160);

  DUT* m164 = new DUT("STA-3G", 164,0);
  //m163->lv1Range(3,10);
  //m163->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-163-masks.txt");
  //m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-163-etacalib.txt");
  config.addDut(m164);

}


void may2009DUTs(TbConfig &config){
  //Configure DUTs used in the may2009 test beam
  DUT* m160 = new DUT("PLANAR", 160,0);
  //Range of acceptable lv1 values, keep loose to avoid false positives in noise scan
  m160->lv1Range(2,8);
  m160->addMasks("calibs/may2009-hotpixelfinder-160-masks.txt");
  //Eta corrections, output from getEtaCorr.
  m160->ecorrs.addEcorr("calibs/may2009-boff-a15-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/may2009-boff-a0-1-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/may2009-boff-a0-2-getetacorr-160-etacalib.txt");
  m160->addToTCalib("calibs/totToQ_160.root");
  m160->getToTcalib()->chargeCorrection(8.117/7.0); //wrong C-low
  config.addDut(m160);
  DUT* m164 = new DUT("STA-3G",164,0);
  m164->lv1Range(2,10);
  m164->addMasks("calibs/may2009-hotpixelfinder-164-masks.txt");
  m164->ecorrs.addEcorr("calibs/may2009-boff-a15-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/may2009-boff-a0-1-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/may2009-boff-a0-2-getetacorr-164-etacalib.txt");
  m164->addToTCalib("calibs/totToQ_164.root");
  config.addDut(m164);
  DUT* m168 = new DUT("FBK-3E7",168,0);
  m168->lv1Range(3,11);
  m168->addMasks("calibs/may2009-hotpixelfinder-168-masks.txt");
  m168->ecorrs.addEcorr("calibs/may2009-boff-a15-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/may2009-boff-a0-1-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/may2009-boff-a0-2-getetacorr-168-etacalib.txt");
  m168->addToTCalib("calibs/totToQ_168.root");
  config.addDut(m168);
  DUT* m172 = new DUT("FBK-3EM5",172,0);
  m172->lv1Range(2,10);
  m172->addMasks("calibs/may2009-hotpixelfinder-172-masks.txt");
  m172->ecorrs.addEcorr("calibs/may2009-boff-a15-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/may2009-bon-a0-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/may2009-bon-a15-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/may2009-boff-a0-1-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/may2009-boff-a0-2-getetacorr-172-etacalib.txt");
  m172->addToTCalib("calibs/totToQ_172.root");
  config.addDut(m172);
}

void may2009DUTsOld1(TbConfig &config){
  //Configure DUTs used in the may2009 test beam
  DUT* m160 = new DUT("PLANAR", 160,0);
  m160->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-160-etacalib.txt");
  m160->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-160-etacalib.txt");
  m160->addToTCalib("calibs/totToQ_160.root");
  m160->getToTcalib()->chargeCorrection(8.117/7.0); //wrong C-low
  config.addDut(m160);
  DUT* m164 = new DUT("STA-3G",164,0);
  m164->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-164-etacalib.txt");
  m164->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-164-etacalib.txt");
  m164->addToTCalib("calibs/totToQ_164.root");
  config.addDut(m164);
  DUT* m168 = new DUT("FBK-3E7",168,0);
  m168->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-168-etacalib.txt");
  m168->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-168-etacalib.txt");
  m168->addToTCalib("calibs/totToQ_168.root");
  config.addDut(m168);
  DUT* m172 = new DUT("FBK-3EM5",172,0);
  m172->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a0-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-boff-a15-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a0-getetacorr-172-etacalib.txt");
  m172->ecorrs.addEcorr("calibs/unprocessed/may2009-unprocessed-bon-a15-getetacorr-172-etacalib.txt");
  m172->addToTCalib("calibs/totToQ_172.root");
  config.addDut(m172);
}

//Configuration for BAT May2009 
void BATmay2009(TbConfig &config){
  cout << "Using configuration may2009" << endl;
  //Site configuration
  #ifdef SITEBATMAY2009_SET
  siteBATMay2009( config );
  #else 
  cout << "Site config for BATmay2009 not found" << endl;
  #endif
  // Set up and configure eventbuilders
  config.addBuilder(new BatTrack);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  Chi2Builder* chi2 = new Chi2Builder(15.0); //Chi2 cut
  config.addBuilder( (EventBuilder*) chi2);
  AngleCuts* angles = new AngleCuts( "calibs/may2009-angledist.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
  config.eventbuilders.push_back( (EventBuilder*) angles );
  DutSync* dutsync = new DutSync();
  dutsync->readFile("calibs/dutsync-checkdutsync-160-dutsync.txt");
  config.addBuilder( dutsync );
  may2009DUTs(config);
}

//Configuration for BAT May2009 tbtrack.1 files
void BATmay2009old1(TbConfig &config){
  cout << "Using configuration may2009.tbtrack.1" << endl;
  //Site configuration
  #ifdef SITEBATMAY2009OLD1_SET
  siteBATMay2009old1( config );
  #else 
  cout << "Site config for BATmay2009old1 not found" << endl;
  #endif
  // Set up and configure eventbuilders
  config.addBuilder(new BatTrack);
  config.addBuilder(new PixelMasker); // need this to avoid empty hits vector!
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  Chi2Builder* chi2 = new Chi2Builder(20.0); //Chi2 cut
  config.addBuilder( (EventBuilder*) chi2);
  AngleCuts* angles = new AngleCuts( "calibs/unprocessed/may2009-unprocessed-angledist.txt", 1.5); //Angle cut at mean +/- 1.5 sigma
  config.addBuilder( (EventBuilder*) angles );
  DutSync* dutsync = new DutSync();
  dutsync->readFile("calibs/unprocessed/run500-1500-checkdutsync-160-dutsync.txt");
  config.addBuilder( dutsync );
  Translator* translator = new Translator();
  //160
  translator->addTranslation(160, 600, 803, 6.0, 1.23e-1);
  translator->addTranslation(160, 804, 998, 12.7, -1.69);
  translator->addTranslation(160, 999,1204, 6.85,5.03e-1);
  translator->addTranslation(160,1205,1307, -3.01, -5.55e-2);
  translator->addTranslation(160,1308,1364, -1.79, 2.11e-1);
  //164
  translator->addTranslation(164, 600, 803, -6.9, 3.36e-1);
  translator->addTranslation(164, 804, 998, 13.1, -1.32e-1);
  translator->addTranslation(164, 999,1204, 3.51, 7.34e-4);
  translator->addTranslation(164,1205,1307, -6.37, 1.93e-1);
  translator->addTranslation(164,1308,1364, -6.81, 1.23e-1);
  //168
  translator->addTranslation(168, 600, 803, -9.58, 2.69e-1);
  translator->addTranslation(168, 804, 998, 22.5, 2.0e-1);
  translator->addTranslation(168, 999,1204, 10.9, -1.38e-2);
  translator->addTranslation(168,1205,1307, -7.38, 3.64e-1);
  translator->addTranslation(168,1308,1364, -7.48, -3.79e-1);
  //172
  translator->addTranslation(172, 600, 803, -5.22, 6.23e-1);
  translator->addTranslation(172, 804, 998, 32.9, -1.534);
  translator->addTranslation(172, 999,1204, -1.40e-1,-1.33);
  translator->addTranslation(172,1205,1307, -1.15, -1.22);
  translator->addTranslation(172,1308,1364, -1.09, -2.75e-1);
  config.addBuilder(translator);
  
  may2009DUTsOld1(config);
}


void BAToct2009(TbConfig &config){
  cout << "Using configuration oct2009" << endl;
  #ifdef SITEBATOCT2009_SET
  siteBATOct2009( config );
  #else 
  cout << "Site config for BAToct2009 not found" << endl;
  #endif
  // Set up and configure eventbuilders
  config.addBuilder(new BatTrack);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  Chi2Builder* chi2 = new Chi2Builder(20.0); //Chi2 cut at 20.
  config.addBuilder( (EventBuilder*) chi2);
  AngleCuts* angles = new AngleCuts( "calibs/bat-oct2009-trackangle.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
  config.addBuilder( (EventBuilder*) angles );
  DutSync* dutsync = new DutSync();
  dutsync->readFile("calibs/Oct09_h8_Bat_all-checkdutsync-163-dutsync.txt");
  config.addBuilder( dutsync );

  DUT* m162 = new DUT("PLANAR", 162,0);
  m162->lv1Range(3,10);
  m162->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-162-masks.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus7degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus15degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus22degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus22degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus30degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus30degree-getetacorr-162-etacalib.txt");
  m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus7degree-getetacorr-162-etacalib.txt");
  config.addDut(m162);

  DUT* m163 = new DUT("STA", 163,0);
  m163->lv1Range(3,10);
  m163->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-163-masks.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus7degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus15degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus22degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus22degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus30degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus30degree-getetacorr-163-etacalib.txt");
  m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus7degree-getetacorr-163-etacalib.txt");
  config.addDut(m163);
  DUT* m165 = new DUT("FBK", 165,0);
  m165->lv1Range(2,10);
  m165->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-165-masks.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus7degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus15degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus22degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus22degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus30degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.minus30degree-getetacorr-165-etacalib.txt");
  m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.plus7degree-getetacorr-165-etacalib.txt");
  config.addDut(m165);
}

//Configuration for BAT Oct 2010
void BAToct2010(TbConfig &config){
  cout << "Using configuration Oct 2010" << endl;
  //Site configuration
  #ifdef SITEBATOCT2010_SET
  siteBATOct2010( config );
  #else 
  cout << "Site config for BAT  Oct 2010 not found" << endl;
  #endif
  // Set up and configure eventbuilders
  config.addBuilder(new BatTrack);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  Chi2Builder* chi2 = new Chi2Builder(10.0); //Chi2 cut
  config.addBuilder( (EventBuilder*) chi2);
  //AngleCuts* angles = new AngleCuts( "calibs/Oct2010-BAT-fieldOFF-0degree-trackangle.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
  AngleCuts* angles = new AngleCuts( "calibs/Oct2010-BAT-trackangle.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
  config.eventbuilders.push_back( (EventBuilder*) angles );
  //DutSync* dutsync = new DutSync();
  //dutsync->readFile("calibs/dutsync-checkdutsync-160-dutsync.txt");
  //config.addBuilder( dutsync );

  //Configure DUTs used in the Oct 2009 test beam
  DUT* m160 = new DUT("STA-3G", 160,0);
  //Range of acceptable lv1 values, keep loose to avoid false positives in noise scan
  m160->lv1Range(6,10);
  m160->addMasks("calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-160-masks.txt");
  //Eta corrections, output from getEtaCorr.
  //m160->addToTCalib("calibs/totToQ_160.root");
  //m160->getTotcalib()->chargeCorrection(8.117/7.0); //wrong C-low
  config.addDut(m160);

  DUT* m161 = new DUT("FBK-3E-5E15",161,0);
  m161->lv1Range(6,9);
  m161->addMasks("calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-161-masks.txt");
  //m161->addToTCalib("calibs/totToQ_164.root");
  config.addDut(m161);

  DUT* m162 = new DUT("FBK-4E-3E15",162,0);
  m162->lv1Range(6,10);
  m162->addMasks("calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-162-masks.txt");
  //m162->addToTCalib("calibs/totToQ_168.root");
  config.addDut(m162);

  DUT* m163 = new DUT("PLANAR",163,0);
  m163->lv1Range(6,9);
  m163->addMasks("calibs/Oct2010-BAT-fieldOFF-0degree-hotpixelfinder-163-masks.txt");
  //m163->addToTCalib("calibs/totToQ_172.root");
  config.addDut(m163);
}

void EudetOct2009(TbConfig& config){
#ifdef SITEEUDETOCT2009_SET
  siteEudetOct2009(config);
#else 
  cout << "Site config for EudetOct2009 not found" << endl;
#endif
  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  config.addBuilder(new CheckRegion);
  eut->addMatchDUT(10);
  eut->addMatchDUT(11);
  eut->addMatchDUT(13);
  eut->nMatches(1);
  void addTranslation(int iden, int firstRun, int lastRun, double shiftX, double shiftY);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  DUT* m10 = new DUT("FBK", 10, 3);
  config.addDut(m10);
  DUT* m11 = new DUT("STA", 11, 3);
  config.addDut(m11);
  DUT* m13 = new DUT("PLA", 13, 3);
  config.addDut(m13);
  DUT* m1 = new DUT("TEL", 1, 0);
  config.addDut(m1);
  DUT* m2 = new DUT("TEL", 2, 0);
  config.addDut(m2);
  DUT* m3 = new DUT("TEL", 3, 0);
  config.addDut(m3);
  DUT* m4 = new DUT("TEL", 4, 0);
  config.addDut(m4);
  DUT* m5 = new DUT("TEL", 5, 0);
  config.addDut(m5);
} 

void EudetPPS2009(TbConfig& config){
#ifdef SITEEUDETOCT2009_SET
  siteEudetOct2009(config);
#else 
  cout << "Site config for Eudet2009 not found" << endl;
#endif
  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  eut->addMatchDUT(10);
  eut->addMatchDUT(14);
  eut->nMatches(1);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  DUT* m10 = new DUT("SIN", 10, 3);
  m10->lv1Range(0,10);
  config.addDut(m10);
  DUT* m15 = new DUT("SIN", 15, 3);
  m15->lv1Range(0,10);
  config.addDut(m15);
  DUT* m14 = new DUT("SIN", 14, 3);
  m14->lv1Range(0,10);
  config.addDut(m14);
}

void EudetNov2009(TbConfig& config){
#ifdef SITEEUDETNOV2009_SET
  siteEudetNov2009(config);
#else 
  cout << "Site config for EudetNov2009 not found" << endl;
#endif
  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  eut->addMatchDUT(14);
  eut->addMatchDUT(15);
  eut->addMatchDUT(17);
  eut->nMatches(0);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  DUT* m14 = new DUT("PLA", 14, 3);
  m14->lv1Range(0,10);
  config.addDut(m14);
  DUT* m15 = new DUT("STA", 15, 3);
  m15->lv1Range(0,10);
  config.addDut(m15);
  DUT* m17 = new DUT("SIN", 17, 3);
  m17->lv1Range(0,10);
  config.addDut(m17);
}


void EudetJune2010(TbConfig &config){
#ifdef SITEEUDETJUNE2010_SET
  siteEudetJune2010(config);
#else 
  cout << "Site config for EudetJune2010 not found" << endl;
#endif
  cout << config.dataPath << endl;
  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  eut->addMatchDUT(10);
  eut->addMatchDUT(11);
  eut->addMatchDUT(12);
  //eut->addMatchDUT(13);
  eut->nMatches(1);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  DUT* m10 = new DUT("PLA", 10, 3);
  m10->lv1Range(1,9);
  m10->addMasks("calibs/june2010-eudet-all-hotpixelfinder-10-masks.txt");
  m10->ecorrs.addEcorr("calibs/june2010-eudet-0degree-getetacorr-10-etacalib.txt");
  m10->ecorrs.addEcorr("calibs/june2010-eudet-15degree-getetacorr-10-etacalib.txt");
  config.addDut(m10);

  DUT* m11 = new DUT("FBK", 11, 3);
  m11->lv1Range(0,10);
  m11->addMasks("calibs/june2010-eudet-all-hotpixelfinder-11-masks.txt");
  m11->ecorrs.addEcorr("calibs/june2010-eudet-0degree-getetacorr-11-etacalib.txt");
  m11->ecorrs.addEcorr("calibs/june2010-eudet-15degree-getetacorr-11-etacalib.txt");
  config.addDut(m11);

  DUT* m12 = new DUT("FBK", 12, 3);
  m12->lv1Range(0,10);
  m12->addMasks("calibs/june2010-eudet-all-hotpixelfinder-12-masks.txt");
  m12->ecorrs.addEcorr("calibs/june2010-eudet-0degree-getetacorr-12-etacalib.txt");
  m12->ecorrs.addEcorr("calibs/june2010-eudet-15degree-getetacorr-12-etacalib.txt");
  config.addDut(m12);

  DUT* m13 = new DUT("FBK", 13, 3);
  m13->lv1Range(2,14);
  //m13->addMasks("calibs/dut13-masks.txt");
  //m13->ecorrs.addEcorr("calibs/june2010_eudet_0degree-getetacorr-13-etacalib.txt");
  //m13->ecorrs.addEcorr("calibs/june2010_eudet_15degree-getetacorr-13-etacalib.txt");
  config.addDut(m13);
}


void EudetNov2010(TbConfig &config){
#ifdef SITEEUDETNOV2010_SET
  siteEudetNov2010(config);
#else 
  cout << "Site config for EudetNov2010 not found" << endl;
#endif
  cout << config.dataPath << endl;
  sprintf(config.treeName, "eutracks");
  EuBuildTrack* eut = new EuBuildTrack();
  //eut->addMatchDUT(10);
  //eut->addMatchDUT(11);
  //eut->addMatchDUT(12);
  //eut->addMatchDUT(13);
  //eut->nMatches(1);
  config.addBuilder(eut);
  config.addBuilder(new PixelMasker);
  config.addBuilder(new ClusterFinder);
  config.addBuilder(new CheckRegion);
  DUT* m10 = new DUT("STA3D", 10, 3);
  m10->lv1Range(1,10);
  //m10->addMasks("calibs/*.txt");
  //m10->ecorrs.addEcorr("calibs/*.txt");
  config.addDut(m10);
  DUT* m11 = new DUT("PPS1", 11, 3);
  m11->lv1Range(1,10);
  //m11->addMasks("calibs/*.txt");
  //m11->ecorrs.addEcorr("calibs/*.txt");
  config.addDut(m11);
  DUT* m12 = new DUT("PPS0", 12, 3);
  m12->lv1Range(1,10);
  //m12->addMasks("calibs/*.txt");
  //m12->ecorrs.addEcorr("calibs/*.txt");
  config.addDut(m12);
}


void SIMbat2010(TbConfig &config){
  cout << "Using configuration SimBat2010" << endl;
  #ifdef SITESIMBAT2010_SET
  siteSimBat2010( config );
  #else 
  cout << "Site config for simBat2010 not found" << endl;
  #endif
  // Set up and configure eventbuilders
  /*config.eventbuilders.push_back(new BatTrack);
  config.eventbuilders.push_back(new PixelMasker);
  config.eventbuilders.push_back(new ClusterFinder());
  Chi2Builder* chi2 = new Chi2Builder(15.0); //Chi2 cut at 20.
  config.eventbuilders.push_back( (EventBuilder*) chi2);
  //AngleCuts* angles = new AngleCuts( "calibs/bat-oct2009-trackangle.txt" , 1.5); //Angle cut at mean +/- 1.5 sigma
  //config.eventbuilders.push_back( (EventBuilder*) angles );
  //DutSync* dutsync = new DutSync();
  //dutsync->readFile("calibs/Oct09_h8_Bat_all-checkdutsync-163-dutsync.txt");
  //config.eventbuilders.push_back( dutsync );

  config.eventbuilders.push_back(new simBaseBuilder());
  //config.eventbuilders.push_back(new simTruthBuilder());
   */

  //Data getters
  config.addBuilder(new BatTrack());
  config.addBuilder(new simBaseBuilder());
  config.addBuilder(new simTruthBuilder());
  config.addBuilder(new simPixelEdepBuilder());
  //Simulator
  config.addBuilder(new simDutRunner());
  //Masks & cuts
  config.addBuilder(new PixelMasker());
  config.addBuilder(new ClusterFinder());
  config.addBuilder(new CheckRegion());
  config.addBuilder(new Chi2Builder(15.0));
  //config.addBuilder(new etaCutter(0.2,0.8));

  DUT* m160 = new DUT("PLANAR", 160,0);
  m160->lv1Range(0,1);
  //m162->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-162-masks.txt");
  //m162->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-162-etacalib.txt");
  config.addDut(m160);
  config.simIdenNameMap[160] = "DUT1";

  DUT* m164 = new DUT("STA-3G", 164,0);
  m164->lv1Range(0,1);
  //m163->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-163-masks.txt");
  //m163->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-163-etacalib.txt");
  config.addDut(m164);
  config.simIdenNameMap[164] = "DUT2";

  DUT* m168 = new DUT("FBK-3E7", 168,0);
  m168->lv1Range(0,1);
  //m165->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-165-masks.txt");
  //m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-165-etacalib.txt");
  config.addDut(m168);
  config.simIdenNameMap[168] = "DUT3";

  DUT* m172 = new DUT("FBK-3EM5", 172,0);
  m172->lv1Range(0,1);
  //m165->addMasks("calibs/Oct09_h8_Bat_all-hotpixelfinder-165-masks.txt");
  //m165->ecorrs.addEcorr("calibs/Oct09_H8_Bat_FieldON.0degree-getetacorr-165-etacalib.txt");
  config.addDut(m172);
  config.simIdenNameMap[172] = "DUT4";

  //Needed for simDataKeeper::getHitsByIden()
  config.simIdenNameMap_all[81] = "BatMod1";
  config.simIdenNameMap_all[83] = "BatMod3";
  config.simIdenNameMap_all[86] = "BatMod6";
  //Copy DUT idens into full list
  for (map<int,string>::iterator it = config.simIdenNameMap.begin();
      it != config.simIdenNameMap.end(); it++) {
    config.simIdenNameMap_all[it->first] = it->second;
  }
}
