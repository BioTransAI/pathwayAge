Project title: PathwayAge

Project description: PathwayAge is a 2-stage biological age predictor that first builds separate machine learning models for methylation sites mapping to each of the BP pathways, yielding 1 machine-learning model per pathway (first-stage). This procedure compresses data from individual methylation sites into a pathway-level feature. Then, a second-stage algorithm integrates these pathway-level features into a systems-level regression.

Real datasets: The datasets in ./PathwayAge/data/ folder are real datasets used in our experiments. 

To explore the pathwayAge model just run 15.phase3.py. To test the code, only 5 pathways were used and the test results in the ./PathwayAge/result/ folder. Please check the path in the yaml file before. 

Any question please contact: 

Have Fun!