$(document).ready(function () {
    // Button click event
    $("button").click(function () {
        $(".overlay").addClass("visible");
        $(".main-popup").addClass("visible");
    });

    // Close popup event
    $("#popup-close-button a, .overlay").click(function () {
        $(".overlay").removeClass("visible");
        $(".main-popup").removeClass("visible");
    });
});
