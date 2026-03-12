HOWTO: Edit and regenerate the align2graph workflow figure
===========================================================

Files
-----
  align2graph_workflow.mmd   Mermaid source diagram (edit this)
  align2graph_workflow.png   Rendered output (publication-quality PNG)


Dependencies
------------
Mermaid CLI (mmdc) is required. Install once via npm:

  npm install -g @mermaid-js/mermaid-cli

Puppeteer (headless Chrome) needs a config file to disable the sandbox
when running as root or in certain Linux environments. Create it once:

  echo '{"args":["--no-sandbox","--disable-setuid-sandbox"]}' > /tmp/mmdc_puppeteer.json


Editing the diagram
-------------------
Open align2graph_workflow.mmd in any text editor.
The file uses Mermaid flowchart syntax (https://mermaid.js.org/syntax/flowchart.html).
Node labels support basic HTML tags: <br/>, <i>, <b>.
Colours and layout are controlled in the %%{init: ...}%% header block at the top.


Regenerating the figure
-----------------------
Run from the miscellaneous/ directory:

  mmdc -i align2graph_workflow.mmd \
       -o align2graph_workflow.png \
       -w 2800 --scale 2 -b white \
       -p /tmp/mmdc_puppeteer.json

Options:
  -w 2800      canvas width in pixels (wide enough to keep long TSV example on one line)
  --scale 2    pixel density multiplier (~5600px effective width; suitable for print)
  -b white     white background

To also generate an SVG (e.g. for web use):

  mmdc -i align2graph_workflow.mmd \
       -o align2graph_workflow.svg \
       -w 2800 -b white \
       -p /tmp/mmdc_puppeteer.json
