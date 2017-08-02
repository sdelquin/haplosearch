$ ->
    $("#process").on("submit", handle_process)
    $("#download").on("submit", handle_download)
    $(window).on("scroll", handle_scroll)
    $("#back-to-top").on("click", back_to_top)
    $("#back-to-top").tooltip("show")


handle_process = (event) ->
    event.preventDefault()
    $("#process_submit_btn").html("
        <span class='btn-label'>
            <i class='fa fa-spinner fa-spin'></i>
        </span>
        Processing your file...
    ")
    this.submit()

handle_download = (event) ->
    event.preventDefault()
    $("#download_submit_btn").prop("disabled", true)
    this.submit()

handle_scroll = ->
    if $(this).scrollTop() > 50
        $("#back-to-top").fadeIn()
    else
        $("#back-to-top").fadeOut()

back_to_top = ->
    $(this).tooltip("hide")
    $(this).animate
        scrollTop: 0
        800
