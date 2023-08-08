//CoI3 filter runing for 1024 cycles
`timescale 1ns/1ns

module COI3 (input wire               clk,
			 input wire               rst_in,
			 input wire        [10:0] N_in,
			 input wire         	  d_in,
			 output wire       [27:0] d_out,
			 output reg               done
			);

// Integrator stage registers
reg [10:0] d1;
reg [19:0] d2;
reg [27:0] d3;
reg [10:0] count;

//======== Synchronize reset for AD ==========
wire rstb_raw    = !rst_in;
reg  rstb_syn1, rstb;
wire sclk_sync  = clk  && !rstb;
always @(negedge sclk_sync or negedge rstb_raw)
begin
    if (!rstb_raw) 
    	begin
        	rstb_syn1 <= 1'b0;
        	rstb      <= 1'b0;
    	end
    else 
    	begin
        	rstb_syn1 <= 1'b1;
        	rstb      <= rstb_syn1;
    	end
end
//======== Synchronize reset for AD ==========

wire gclk = clk && !done && rstb;
always @(negedge gclk or negedge rstb)
begin
	if (!rstb)
		begin
			d1 <= 0;
			d2 <= 0;
			d3 <= 0;				
			count  <= 0;
			done   <= 0;
		end 
	else
		begin
		// Integrator section
			d1 <= d_in + d1;				
			d2 <= d1   + d2;				
			d3 <= d3   + d2;							
			// Decimation
			if (count == (N_in+1))  //one cycle delay	
				done <= 1'b1;				
			count <= count + 1'b1;
		end
end

assign d_out = done ? d3 : 0;					
endmodule
