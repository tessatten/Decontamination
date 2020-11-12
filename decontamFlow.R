

library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)

decontamFlow <- DiagrammeR::grViz("digraph {

  graph [layout = dot, color = black, fillcolor = Gainsboro, rankdir = TB, style = 'filled, rounded, bold']
  
  node[style='filled, bold, rounded', shape = box, fillcolor = CornSilk, fontsize = 30]
  
  edge[arrowhead = vee, color = black] 

  # nodes
  start [label='All genera', group = g1]
  notContam [label='NOT indicated as likely \ncontaminants using\n Decontam or positive controls: \nretained', fillcolor = Honeydew]
  firstContam [label='Indicated as likely \ncontaminants using \nDecontam', group = g1]
  posContam [label='Indicated as likely \ncontaminants using \n positive controls']
  allContam [label='All likely contaminants', group = g1]
  inGut [label='Previously identified\n in the human \ngut/stool: retained', fillcolor = Honeydew]
  notGut [label='NOT previously identified\n in the human \ngut/stool', group = g1]
  contamorNoHuman [label='Reported as a \ncontaminant OR no \nknown link to \nhumans', group = g1]
  ecologyNoHuman [label='Unlikely ecology\n AND no known\n link to humans: \nremoved', fillcolor = MistyRose]
  contamNoHuman [label='Reported as a \ncontaminant AND no \nknown link to \nhumans: removed', fillcolor = MistyRose]
  noNegs [label='In no negative\n extractions: retained', fillcolor = Honeydew]
  medNegs [label='In one or two\n negative extractions', group = g1]
  manyNegs [label='In three or more\n negative extractions: \nremoved', fillcolor = MistyRose]
  ecologyLikely [label='Ecology indicates plausible\n presence in human gut: \nretained', fillcolor = Honeydew]
  ecologyunlikely [label='Ecology indicates unlikely\n in human gut: \nremoved', fillcolor = MistyRose]


  
  # order
  start -> notContam
  start -> posContam
  start -> firstContam 
  posContam -> allContam
  firstContam -> allContam
  allContam -> inGut
  allContam -> notGut
  notGut -> contamNoHuman
  notGut -> contamorNoHuman
  notGut -> ecologyNoHuman
  contamorNoHuman -> noNegs
  contamorNoHuman -> manyNegs
  contamorNoHuman -> medNegs
  medNegs -> ecologyLikely
  medNegs ->ecologyunlikely
  
  subgraph cluster0 {
  graph [fontsize = 32]
  edge[arrowhead = none, color = Gainsboro] 
  label = 'Key';
  removed [label = 'Removed', fillcolor = MistyRose, group = g2]
  process [label = 'Undetermined', fillcolor = CornSilk, group = g2]
  retained [label = 'Retained', fillcolor = Honeydew, group = g2]


  }
}")

decontamFlow

#ranksep = 0.3, nodesep = 0.3, 


#####add in the names of genera

namesdecontamFlow <- DiagrammeR::grViz("digraph {

  graph [layout = dot, color = black, fillcolor = Gainsboro, rankdir = TB, style = 'filled, rounded, bold']
  
  node[style='filled, bold, rounded', shape = box, fillcolor = CornSilk, fontsize = 30]
  
  edge[arrowhead = vee, color = black] 

  # nodes
  start [label='All genera (195)', group = g1]
  notContam [label='NOT indicated as likely \ncontaminants using\n Decontam or positive controls: \nretained (149)', fillcolor = Honeydew]
  firstContam [label='Indicated as likely \ncontaminants using \nDecontam (46)', group = g1]
  posContam [label='Indicated as likely \ncontaminants using \n positive controls']
  allContam [label='All likely contaminants', group = g1]
  inGut [label='Previously identified\n in the human \ngut/stool: retained (17)', fillcolor = Honeydew]
  notGut [label='NOT previously identified\n in the human \ngut/stool (29)', group = g1]
  contamorNoHuman [label='Reported as a \ncontaminant OR no \nknown link to \nhumans (13)', group = g1]
  ecologyNoHuman [label='Unlikely ecology\n AND no known\n link to humans: \nremoved (4)', fillcolor = MistyRose]
  contamNoHuman [label='Reported as a \ncontaminant AND no \nknown link to \nhumans: removed (12)', fillcolor = MistyRose]
  noNegs [label='In no negative\n extractions: retained (3)', fillcolor = Honeydew]
  medNegs [label='In one or two\n negative extractions (8)', group = g1]
  manyNegs [label='In three or more\n negative extractions: \nremoved (2)', fillcolor = MistyRose]
  ecologyLikely [label='Ecology indicates plausible\n presence in human gut: \nretained (5)', fillcolor = Honeydew]
  ecologyunlikely [label='Ecology indicates unlikely\n in human gut: \nremoved (3)', fillcolor = MistyRose]


  
  # order
  start -> notContam
  start -> posContam
  start -> firstContam 
  posContam -> allContam
  firstContam -> allContam
  allContam -> inGut
  allContam -> notGut
  notGut -> contamNoHuman
  notGut -> contamorNoHuman
  notGut -> ecologyNoHuman
  contamorNoHuman -> noNegs
  contamorNoHuman -> manyNegs
  contamorNoHuman -> medNegs
  medNegs -> ecologyLikely
  medNegs ->ecologyunlikely
  
  subgraph cluster0 {
  graph [fontsize = 32]
  edge[arrowhead = none, color = Gainsboro] 
  label = 'Key';
  removed [label = 'Removed', fillcolor = MistyRose, group = g2]
  process [label = 'Undetermined', fillcolor = CornSilk, group = g2]
  retained [label = 'Retained', fillcolor = Honeydew, group = g2]


  }
}")

namesdecontamFlow