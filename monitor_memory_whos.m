function [ memory_in_use ] = monitor_memory_whos( )
%MONITOR_MEMORY_WHOS uses the WHOS command and evaluates inside the BASE
%workspace and sums up the bytes.  The output is displayed in MB.

mem_elements = evalin('base','whos');
if size(mem_elements,1) > 0

    for i = 1:size(mem_elements,1)
        memory_array(i) = mem_elements(i).bytes;
    end

    memory_in_use = sum(memory_array);
    memory_in_use = memory_in_use/1048576;
else
    memory_in_use = 0;
end
