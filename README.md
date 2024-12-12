main is written in BiCNDP.py (has the main model and constraints)
Callbacks are implemented in BiCNDP_callbacks.py
Instance generator class (generator and read from and write into file methods) are in geninstance.py
Requires a parameter parfile.txt 
  number of nodes, graph type (er or ws), p is ER model probability parameter and if WS it's the rewiring probability parameter 
  edge_density is K parameter of WS model
  ws_param is a parameter for rewiring std dev (rewiring is randomly chosen ~ N(p, ws_param)
  bugetprop is the proportion of nodes chosen as budget (same for both leader and attacker atm - should be changed)
  dbug mode with more verbosity
  test mode runs on a toy wheel graph
  ccut toggles connectivity cuts on and off
  os mac (m) or linux
  file 1 - for writing the generated instance into the data folder (without running). 2- reads the file and runs the model
