getMicroSyntenySearch <- function() {
  paste(
    "$.ajax({",
      "url: 'https://' + url2 + '/services/v1/micro-synteny-search/',",
      "dataType: 'json',",
      "data: JSON.stringify({",
        "query: families1,",
        "matched: $('input#matched').val(),",
        "intermediate: $('input#intermediate').val()",
      "}),",
      "type: 'POST',",
      "success: function(response2) {",
        "obj2 = response2;",
        # Send information about neighboring and related genes back to the chart
        "Shiny.onInputChange('genomicLinkages', {",
          "results1: obj1,",
          "results2: obj2",
        "});",
      "},",
      "error: function(errmsg2) { alert('FAIL2: ' + errmsg2.responseText); }",
    "});"
  )
}

getGeneToQueryTrack <- function(url1, url2, microSyntenySearch) {
  paste(
    # Query the Genome Context Viewer for genomic linkages
    sprintf("url1 = '%s';", url1),
    sprintf("url2 = '%s';", url2),
    "var geneString = '';",
    "mt0 = this.url.search('medtr');",
    "cc0 = this.url.search('cajca');",
    "gm0 = this.url.search('glyma');",
    "vu0 = this.url.search('vigun');",
    "js1 = this.url.search('/json');",
    "if (mt0 >= 0) {",
      "geneString = this.url.substring(mt0, js1);",
    "} else if (gm0 >= 0) {",
      "geneString = this.url.substring(gm0, js1);",
    "} else if (vu0 >= 0) {",
      "geneString = this.url.substring(vu0, js1);",
    "} else if (cc0 >= 0) {",
      "geneString = this.url.substring(cc0, js1);",
    "} else if (this.url.search('arabidopsis.org') >= 0) {",
      "at0 = this.url.search('name=');",
      "geneString = 'arath.Col.' + this.url.substring(at0 + 5);",
    "} else {",
      "return;",
    "}",
    "$.ajax({",
      "url: 'https://' + url1 + '/services/v1/gene-to-query-track/',",
      "dataType: 'json',",
      "data: JSON.stringify({",
        "gene: geneString,",
        "neighbors: $('input#neighbors').val()",
      "}),",
      "type: 'POST',",
      "success: function(response) {",
        "obj1 = response;",
        "families1 = Array.from(obj1.genes, x => x.family);",
        microSyntenySearch,
      "},",
      "error: function(errmsg) { alert('FAIL: ' + errmsg.responseText); }",
    "});"
  )
}

getProvideMultipleURLs <- function() {
  paste(
    # From the JSON at this.url, extract the URLs related to this gene.
    # Note that this.url = legumeInfo_urlBase + geneString + '/json'
    #  legumeInfo_urlBase currently has 34 characters (see above)
    #  and geneString = <5-character species abbreviation>.geneName
    # And for now, add the gene family phylogram URL by hand.
    "$.getJSON(this.url, function(data) {",
      "var geneString = this.url.substring(34, this.url.indexOf('/json'));",
      "var geneName = geneString.substring(6);",
      "var content = '';",
      "if (data.length == 0) {",
        "content = '<p>No ' + geneName + ' links found.</p>';",
      "} else {",
        "$.each(data, function(i, obj) {",
          "content = content + '<p><a href=' + obj.href + ' target=_blank>' + obj.text + '</a></p>';",
          "if (i == 0) {",
            "var urlPhylogram = 'http://legumeinfo.org/chado_gene_phylotree_v2?gene_name=' + geneString;",
            "var textPhylogram = 'View LIS gene family phylogram page for : ' + geneName;",
            "content = content + '<p><a href=' + urlPhylogram + ' target=_blank>' + textPhylogram + '</a></p>';",
          "}",
        "});",
      "}",
      "var $div = $('<div></div>');",
      "$div.html(content);",
      "$div.dialog({",
        "title: geneName + ' Links',",
        "width: 512,",
        "height: 'auto',",
        "modal: true",
      "});",
    "});"
  )
}
