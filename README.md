**FILES AND PARAMETERS**
<p>main is written in BiCNDP.py (has the main model and constraints)</p>
<p>Callbacks are implemented in BiCNDP_callbacks.py</p>
<p>Instance generator class (generator and read from and write into file methods) are in geninstance.py</p>
<p>Requires a parameter parfile.txt </p>
<ul>
  <li>number of nodes, graph type (er or ws)</li>
  <li>p is ER model probability parameter and if WS it's the rewiring probability parameter </li>
  <li>edge_density is K parameter of WS model</li>
  <li>ws_param is a parameter for rewiring std dev (rewiring is randomly chosen ~ N(p, ws_param)</li>
  <li>bugetprop is the proportion of nodes chosen as budget (same for both leader and attacker atm - should be changed)</li>
  <li>dbug mode with more verbosity</li>
  <li>test mode runs on a toy wheel graph</li>
  <li>ccut toggles connectivity cuts on and off</li>
  <li>os mac (m) or linux</li>
  <li>file 1 - for writing the generated instance into the data folder (without running). 2- reads the file and runs the model</li>
</ul>
<p>bat.sh is written to create/run (dependingfile parameter 1/2) instances with varying parameters. Inner loop creates 10 random instances for each parameter setting</p>
