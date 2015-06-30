data = {}

$(function () {
    d3.json("/datafiles/test.json", function(error, value){
        console.log(arguments)
        data = value
        layoutProtein(data)
    })
})

var margin = {top: 5, right: 40, bottom: 20, left: 120},
    width = 960 - margin.left - margin.right,
    height = 1880 - margin.top - margin.bottom;


function layoutProtein(data){
    var svg = d3.select("body").append("svg").attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom).append('g')
    
    var graphicContainer = svg.append("g").attr("id", "graphic-container")
    var entered = graphicContainer.selectAll("g").data(data).enter().append(
      "g").append(
      "rect").attr(
      "x", function(d){return d.start_position}).attr('y', function(d, i){
        return 50 + i * 16
      }).attr(
      "width", function(d){ return (d.end_position - d.start_position) * 4 }).attr(
      "height", 16).attr('class', "subsequence").attr('protein', function(d){
        return d.protein_id})
}
