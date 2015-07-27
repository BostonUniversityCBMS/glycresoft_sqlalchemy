var ajaxForm, setupAjaxForm;

ajaxForm = function(formHandle, success, error) {
  $(formHandle).on('submit', function(event) {
    var ajaxParams, data, encoding, handle, method, url;
    event.preventDefault();
    handle = $(this);
    url = handle.attr('action');
    method = handle.attr('method');
    data = new FormData(this);
    encoding = handle.attr('enctype') || 'application/x-www-form-urlencoded; charset=UTF-8';
    ajaxParams = {
      'url': url,
      'method': method,
      'data': data,
      'processData': false,
      'contentType': false,
      'success': success,
      'error': error
    };
    $.ajax(ajaxParams);
  });
};

setupAjaxForm = function(sourceUrl, container) {
  var isModal;
  container = $(container);
  isModal = container.hasClass('modal');
  $.get(sourceUrl).success(function(doc) {
    if (isModal) {
      container.find('.modal-content').html(doc);
      container.openModal();
      container.find('form').submit(function(event) {
        container.closeModal();
      });
    } else {
      container.html(doc);
    }
  });
  container.find('script').each(function(i, tag) {
    var srcURL;
    tag = $(tag);
    srcURL = tag.attr('src');
    if (srcURL !== void 0) {
      $.getScript(srcURL);
    } else {
      eval(tag.text());
    }
  });
};

//# sourceMappingURL=ajax-form.js.map
