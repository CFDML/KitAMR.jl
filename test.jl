# # input_channel = Channel{String}(1)

# # # 异步任务，用于监听键盘输入
# # function listen_for_input()
# #     @async begin
# #         while true
# #             input = readline(stdin)
# #             put!(input_channel, input)
# #         end
# #     end
# # end

# # # 主任务，用于执行主要逻辑并查询输入
# # function main_task()
# #     i = 0
# #     while true
# #         println("正在执行主任务: $i")
        
# #         # 检查输入通道是否有新的输入
# #         if isready(input_channel)
# #             input = take!(input_channel)
# #             println("收到输入：$input")
# #             # 根据输入执行相应操作
# #             if input=="break"
# #                 break
# #             end
# #         end
# #         i+=1
# #         # sleep(1)  # 模拟其他任务的执行
# #     end
# # end

# # # 启动异步监听键盘输入
# # listen_for_input()

# # # 执行主任务
# # main_task()


# # input_channel = Channel{String}(1)
# # test = Ref(false)
# # function listen_for_input()
# #     @async begin
# #         while true
# #             input = readline(stdin)
# #             put!(input_channel, input)
# #         end
# #     end
# # end
# # function main_task()
# #     i = 0
# #     global test
# #     rank = MPI.Comm_rank(MPI.COMM_WORLD)
# #     while true
# #         input=nothing
# #         if rank==0
# #             @show i
# #             if isready(input_channel)
# #                 input = take!(input_channel)
# #             end
# #             if input=="save"
# #                 test[] = true
# #             end
# #         end
# #         i+=1
# #         MPI.Bcast!(test,0,MPI.COMM_WORLD)
# #         for i = 0:MPI.Comm_size(MPI.COMM_WORLD)-1
# #             if MPI.Comm_rank(MPI.COMM_WORLD)==i
# #                 @show test[]
# #             end
# #             MPI.Barrier(MPI.COMM_WORLD)
# #         end
# #         sleep(1)
# #     end
# # end
# # using MPI
# # MPI.Init()
# # MPI.Comm_rank(MPI.COMM_WORLD)==0 && listen_for_input()
# # main_task()

# using MPI

# input_channel = Channel{String}(1)
# test = Ref(false)

# function listen_for_input()
#     Threads.@spawn begin
#         # println("Rank 0: Listening for input...")
#         while true
#             input = readline(stdin)
#             # println("Rank 0: Read input: $input")
#             put!(input_channel, input)
#         end
#     end
# end

# function main_task()
#     i = 0
#     global test
#     rank = MPI.Comm_rank(MPI.COMM_WORLD)
#     println("Starting main task on rank $rank")

#     while true
#         input = nothing
#         if rank == 0
#             println("Rank 0: Iteration $i")
#             if isready(input_channel)
#                 input = take!(input_channel)
#                 println("Rank 0: Received input: $input")
#             end
#             if input == "save"
#                 test[] = true
#                 println("Rank 0: Setting test to true")
#             end
#         end
#         i += 1
#         MPI.Bcast!(test, 0, MPI.COMM_WORLD)
#         for j = 0:MPI.Comm_size(MPI.COMM_WORLD) - 1
#             if MPI.Comm_rank(MPI.COMM_WORLD) == j
#                 println("Rank $j: test[] = $(test[])")
#             end
#             MPI.Barrier(MPI.COMM_WORLD)
#         end
#         sleep(1)
#     end
# end

# MPI.Init()
# if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#     listen_for_input()
# end
# main_task()

i = 0
while true
    global i
    Base.start_reading(stdin)
    if bytesavailable(stdin)!=0
        input = readline(stdin)
        @show input
    end
    @show i
    i+=1
    sleep(1)
end