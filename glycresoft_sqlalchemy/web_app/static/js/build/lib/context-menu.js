var contextMenu;

contextMenu = function(target, options) {
  var action, item;
  for (item in options) {
    action = options[item];
    console.log(item, action);
  }
  $(target).off("contextmenu", false);
  $(document).off("mousedown", false);
  return $(target).on("contextmenu", function(event) {
    var handle;
    event.preventDefault();
    handle = $(".context-menu");
    handle.empty();
    for (item in options) {
      action = options[item];
      handle.append($("<li></li>").text(item).attr("data-action", item));
    }
    $(".context-menu li").click(function(e) {
      handle = $(this);
      console.log(handle);
      action = options[handle.attr("data-action")];
      return action(handle);
    });
    return $(".context-menu").finish().toggle(100).css({
      top: event.pageY + 'px',
      left: event.pageX + 'px'
    });
  });
};

$(document).on("mousedown", function(e) {
  console.log("Off!");
  if (!$(e.target).parents(".context-menu").length > 0) {
    return $(".context-menu").hide(100);
  }
});

//# sourceMappingURL=context-menu.js.map
