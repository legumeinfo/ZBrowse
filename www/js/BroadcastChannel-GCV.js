// Broadcast Channel allows interapplication communication between web pages (etc)
// with the same origin (protocol, domain, and port).
// https://developer.mozilla.org/en-US/docs/Web/API/Broadcast_Channel_API
// If running ZBrowse locally, note that 127.0.0.1 and localhost are not considered "same origin",
// either will work but they must match for both applications.

// Connect to or disconnect from the LIS Genome Context Viewer
// https://github.com/legumeinfo/lis_context_viewer/wiki/Inter-Application-Communication

// Shiny JavaScript events
// https://shiny.rstudio.com/articles/js-events.html
// Should we connect on shiny:connected or shiny:sessioninitialized?
// It probably does not matter, as messages will be ignored unless a session exists.
$(document).on('shiny:connected', function(e) {
//$(document).on('shiny:sessioninitialized', function(e) {
  try {
    bc = new BroadcastChannel('lis');

    bc.onmessage = function(e) {
      Shiny.onInputChange("bc_gcv", e.data);
      //console.log(e.data);
    };

    //console.log('Opened GCV broadcast channel');
  } catch (ex) {
    //console.log('Could not open a GCV broadcast channel');
  }
});

$(document).on('shiny:disconnected', function(e) {
  bc.close();
  //console.log('Closed GCV broadcast channel');
});
