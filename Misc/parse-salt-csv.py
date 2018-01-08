#!/usr/bin/env python

import csv

#   '#', VEGI ID, Locus, Stock ID, Note, PP1/PP2, Bar rep 1, Bar rep 2, Bar rep 3, Comments,,,,...
def body():
    # <tr><td>2</td><td>AT3G24255</td><td>CS860296</td><td><a href="CS860296" rev="10101008" name="10102008" title="10103008">load</a></td><td><img id="10101008" alt="10101008"></td><td><img id="10102008" alt="10102008"></td><td><img id="10103008" alt="10103008"></td></tr>
    bodyHTML="          <tr><td>{0}</td><td>{1}</td><td>{2}</td><td><a href='{2}' rev='{3}' name='{4}' title='{5}'>load</a></td><td><img id='{3}' alt='{3}'></td><td><img id='{4}' alt='{4}'></td><td><img id='{5}' alt='{5}'></td></tr>"
    with open('Salt_assay_Dec2012.csv', 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        reader.next() # skip header
        for row in reader:
            print bodyHTML.format(row[1], row[2], row[3], row[6], row[7], row[8])
            
            
            
def header():
    print '''
<!DOCTYPE html>
    
<style type="text/css">
    table {
        border-collapse: collapse;
    }
    td, th {
        border: 1px solid black;
        padding: 3px;
    }
    tr {
        text-align: center;
    }
    img {
        width: 700px;
        height: 600px;
    } 
</style>

<html>
    <head>
        <title>Salt Assay Images</title>
        <meta name="author" content="Douglas Hoen">
        <script src="http://code.jquery.com/jquery-latest.js"></script>
        <script type="text/javascript">
            var baseURL = "http://mustang.biol.mcgill.ca:8085/LemnaTec/assay-2012-12/ImagesCollectionVIS/";
            $(document).ready(function() {
                $("a").click(function(event) {
                    event.preventDefault();
                    var barcode1 = this.rev;   // hack
                    var barcode2 = this.name;  // hack
                    var barcode3 = this.title; // hack
                    $('#'+barcode1).html($("<img>").attr("src", baseURL + barcode1 + ".png"));
                    $('#'+barcode2).html($("<img>").attr("src", baseURL + barcode2 + ".png"));
                    $('#'+barcode3).html($("<img>").attr("src", baseURL + barcode3 + ".png"));
                }
            )});
        </script>   
    </head>
    <body>
        <h1>Salt Assay Images</h1>
        Note: Image loading only works in <b>Firefox</b>.<br><br>
        <table>
            <tr><td>VEGI ID</td><td>Locus ID</td><td>Stock ID</td><td>Load Images</td><td>Replicate 1</td><td>Replicate 2</td><td>Replicate 3</td></tr>'''



def footer():
    print '''
        </table>
    </body>
</html>'''



header()
body()
footer()