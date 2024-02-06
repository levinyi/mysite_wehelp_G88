$(document).ready(function() {
    // 获取 CSRF 令牌
    function getCookie(name) {
        let cookieValue = null;
        if (document.cookie && document.cookie !== '') {
            const cookies = document.cookie.split(';');
            for (let i = 0; i < cookies.length; i++) {
                const cookie = cookies[i].trim();
                if (cookie.substring(0, name.length + 1) === (name + '=')) {
                    cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                    break;
                }
            }
        }
        return cookieValue;
    }

    const csrftoken = getCookie('csrftoken');

    // 事件委托，确保动态生成的元素也能绑定事件
    $(document).on('click', '.delete-item-btn', function() {
        var deleteUrl = $(this).data('url');
        var itemId = $(this).data('id'); // 获取项目的ID

        // 显示模态框
        $('#deleteConfirmationModal').modal('show');

        $('#confirmDelete').off('click').on('click', function() {
            $.ajax({
                url: deleteUrl,
                type: "POST",
                data: {
                    'csrfmiddlewaretoken': csrftoken,
                    'id': itemId // 确保也传递项目ID
                },
                dataType: "json",
                success: function (data) {
                    if (data.status === 'success') {
                        alert('Delete successfully!');
                        window.location.reload();
                    } else {
                        alert(data.message);
                    }
                },
                error: function (xhr, errmsg, err) {
                    alert('Delete failed!');
                }
            });

            $('#deleteConfirmationModal').modal('hide');
        });
    });
});
