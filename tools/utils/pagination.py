from django.utils.safestring import mark_safe

class Pagination(object):
    def __init__(self, request, queryset, page_size=10, page_param="page",plus=5):
        """
        :param request: 请求的对象
        :param queryset: 查询的符合条件的数据(根据这个数据给他进行分页处理)
        :param page_size: 每页显示多少条数据
        :param page_param: 在url中传递的获取分页的参数,例如:/number/list/?page=12
        :param plus: 显示当前页的 前或后几页(页码)
        """
        
        import copy
        query_dict = copy.deepcopy(request.GET)
        query_dict._mutable = True
        self.query_dict = query_dict

        self.page_param = page_param
        page = request.GET.get(page_param, "1")
        if page.isdecimal():
            page = int(page)
        else:
            page = 1

        self.page = page
        self.page_size = page_size

        self.start = (page-1) * page_size
        self.end = page * page_size

        self.page_queryset = queryset[self.start: self.end]
        total_count = queryset.count()
        total_page_count, div = divmod(total_count, page_size)
        if div:
            total_page_count += 1
        self.total_page_count = total_page_count
        self.plus = plus

    def html(self):
        if self.total_page_count <= 2*self.plus + 1:
            # 数据库中的数据比较少, 都没有达到11页
            start_page = 1
            end_page = self.total_page_count + 1
        else:
            # 数据库中的数据比较多,>11页.
            # 当前页小于5:
            if self.page <= self.plus:
                start_page = 1
                end_page = 2*self.plus + 1
            else:
                # 当前页>5
                # 当前页+5大于总页码了,那就不合适了
                if (self.page + self.plus) > self.total_page_count:
                    start_page = self.total_page_count - 2*self.plus
                    end_page = self.total_page_count + 1
                else:
                    start_page = self.page - self.plus
                    end_page = self.page + self.plus + 1

        # 页码
        page_str_list = []
        self.query_dict.setlist(self.page_param,[1])
        page_str_list.append('<li class="page-item"><a class="page-link" href="?{}">首页</a></li>'.format(self.query_dict.urlencode()))
        # 上一页
        if self.page >1:
            self.query_dict.setlist(self.page_param,[self.page -1])
            prev = '<li class="page-item"><a class="page-link" href="?{}">上一页</a></li>'.format(self.query_dict.urlencode())
        else:
            self.query_dict.setlist(self.page_param,[1])
            prev = '<li class="page-item"><a class="page-link" href="?{}">上一页</a></li>'.format(self.query_dict.urlencode())
        page_str_list.append(prev)
        
        # 页面
        for i in range(start_page, end_page):
            self.query_dict.setlist(self.page_param,[i])
            if i == self.page:
                ele = '<li class="page-item active"><a class="page-link" href="?{}">{}</a></li>'.format(self.query_dict.urlencode(), i)
            else:
                ele = '<li class="page-item"><a class="page-link" href="?{}">{}</a></li>'.format(self.query_dict.urlencode(), i)
            page_str_list.append(ele)
        
        # 下一页
        if self.page < self.total_page_count:
            self.query_dict.setlist(self.page_param,[self.page + 1])
            next = '<li  class="page-item"><a class="page-link" href="?{}">下一页</a></li>'.format(self.query_dict.urlencode())
        else:
            self.query_dict.setlist(self.page_param,[self.total_page_count])
            next = '<li class="page-item"><a class="page-link" href="?{}">下一页</a></li>'.format(self.query_dict.urlencode())
        page_str_list.append(next)

        # 尾页
        self.query_dict.setlist(self.page_param,[self.total_page_count])
        page_str_list.append('<li class="page-item"><a class="page-link" href="?{}">尾页</a></li>'.format(self.query_dict.urlencode()))
        ### search 
        page_string = mark_safe("".join(page_str_list))
        return page_string