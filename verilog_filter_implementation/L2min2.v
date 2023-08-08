//Function: L2min2 reconstruction filter for IDCs
//No digital trunction applied
//All bitwidths are designed to operate the filter for 1024 cycles (N is 1024)
//The user can change the bitwidth as needed
//Detailed filter information can be found in the paper
//Contact: Dr. Bo Wang, bwang@hbku.edu.qa

`timescale 1ns/1ns

module L2min2 (
			 input wire          		clk,          	//master clock of the filter
			 input wire               	rst_in,	    	//filter reset signal
			 input wire        [10:0]  	N_in, 			//check the paper for the definition of N
			 input wire         	  	d_in,			//1-bit input bitstream
			 output reg       [28:0]  	d_out,			//ceil(3*log2(N)-1) --> 29bit to meet all conditions
			 output wire                done			//reconstruction done signal
			);

// define the counter stage registers
reg [10:0]  d1;  										//bitwidth of ceil(log2(M+1)) --> 11bit to meet all conditions
// define the subtractor stage register
reg [19:0] d2;											//bitwidth of ceil(log2(M^2+M)-1) --> 20bit to meet all conditions

//======== synchronize the reset signal ==========
wire rstb_raw    = !rst_in;
reg  rstb_syn1, rstb;
wire sclk_sync  = clk  && !rstb;
always @(negedge sclk_sync or negedge rstb_raw)
    if (!rstb_raw) begin
        rstb_syn1 <= 1'b0;
        rstb      <= 1'b0;
    end
    else begin
        rstb_syn1 <= 1'b1;
        rstb      <= rstb_syn1;
    end

//======== prepare the weights at clock rising edge ========
wire rst  	  = !rstb;
wire gclk 	  = clk && !done && rstb;
always @(posedge gclk or posedge rst)
begin
	if (rst) 
		d1 <= 0;
	else 
        	d1 <= d1 + 1;
end
//======== the subtractor block operation ========
always @(posedge gclk or posedge rst)
begin
	if (rst)
			d2 <= 20'd524800;  // initial subtractor reset value N_in*(N_in+1)/2, with N_in = 1024 in this example
	else
		begin
			d2 <= d2 - d1;  	//S/A as subtractor
		end
end

//======== last adder stage with clock falling edge control ========
wire clk_add = gclk && d_in;
always @(negedge clk_add or negedge rstb) 
begin
	if (!rstb)
		d_out <= 0;							
	else
		d_out <= d_out + d2;
end

assign done = (d1 == N_in);

endmodule
