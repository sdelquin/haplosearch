$ ->
    $("#process").on("submit", handle_process)
    $("#download").on("submit", handle_download)

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
