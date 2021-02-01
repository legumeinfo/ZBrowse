# Services API v2 calls used in ZZBrowse
# https://github.com/legumeinfo/gcv/wiki/Services-API-v2

microSyntenySearchService <- function(url, families, matched, intermediate) {
  paste(
    sprintf("url = '%s';", url),
    sprintf("families = %s;", families),
    sprintf("matched = %d;", matched),
    sprintf("intermediate = %d;", intermediate),
    "$.ajax({",
      "url: 'https://' + url + '/services/v2/micro-synteny-search/',",
      "dataType: 'json',",
      "data: JSON.stringify({",
        "query: families,",
        "matched: matched,",
        "intermediate: intermediate",
      "}),",
      "type: 'POST',",
      "success: function(response) {",
        "Shiny.onInputChange('microSyntenySearchResults', {",
          "results: response",
        "});",
      "},",
      "error: function(errmsg) {",
        "alert('microSyntenySearchService() failed: ' + errmsg.responseText);",
      "}",
    "});"
  )
}

genesService <- function(url, genes) {
  paste(
    sprintf("genes = %s;", genes),
    sprintf("url = '%s';", url),
    "$.ajax({",
      "url: 'https://' + url + '/services/v2/genes/',",
      "dataType: 'json',",
      "data: JSON.stringify({",
        "genes: genes",
      "}),",
      "type: 'POST',",
      "success: function(response) {",
        "Shiny.onInputChange('genesResults', {",
          "results: response",
        "});",
      "},",
      "error: function(errmsg) {",
        "alert('genesService() failed: ' + errmsg.responseText);",
      "}",
    "});"
  )
}

# Not a Services API call, but convenient to place here as it is
# also a response to a doClickOnLine in the zChart.
provideMultipleURLs <- function(includeGenomicLinkage) {
  paste(
    # From the JSON at this.url, extract the URLs related to this gene.
    # Note that this.url = 'https://legumeinfo.org/gene_links/' + geneString + '/json'
    #  (the first part of which has 34 characters)
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
        # Add the Genomic Linkage button if requested
        ifelse(!includeGenomicLinkage, "", paste(
          "content = content + '<p><button",
            # button actions: first close the dialog, then handle genomic linkages
            "onclick=\"$(this).closest(&quot;.ui-dialog-content&quot;).dialog(&quot;close&quot;);",
            "Shiny.onInputChange(&quot;selectedGene&quot;, &quot;' + geneString + '&quot;);\"",
          ">Genomic Linkage</button></p>';"
        )),
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
