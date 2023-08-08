//Function: L2min2s reconstruction filter for IDCs
//No digital trunction applied
//All bitwidths are designed to operate the filter for maximium 1024 cycles (maximum M is 512)
//The user can change the bitwidth as needed
//Detailed filter information can be found in the paper
//Contact: Dr. Bo Wang, bwang@hbku.edu.qa

`timescale 1ns/1ns

module L2min2s (
			 input wire          		clk,          	//master clock of the filter
			 input wire               	rst_in,	    	//filter reset signal
			 input wire        [9:0]  	M_in, 			//check the paper for the definition of M, in total 2*M-1 sample used for output reconstruction
			 input wire         	  	d_in,			//1-bit input bitstream
			 output reg       [26:0]  	d_out,			//ceil(3*log2(M)) --> 27bit to meet all conditions
			 output wire                done			//reconstruction done signal
			);

// define the counter stage registers
reg [9:0]  d1;  										//bitwidth of ceil(log2(M+1)) --> 10bit to meet all conditions
// define the S/A stage register
reg [17:0] d2;											//bitwidth of ceil(log2(M^2+M)-1) --> 18bit to meet all conditions

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
//reg  state; 					//to indicate the transition from down-counting (0) to up-counting (1) state
reg state;
wire clk_state 	  = (d1 == 0);
always @(posedge clk_state or posedge rst)
begin
	if (rst)
		state 		<= 0;
	else 
		state 		<= 1; 
end	

//always @(posedge gclk or posedge rst)
//begin
//	if (rst)
//		state 		<= 0;
//	else 
//		begin
//		if (d1 == 1) 
//			state 		<= 1; 
//		else 
//			state		 <= state;
//		end
//end	

always @(posedge gclk or posedge rst)
begin
	if (rst) 
		d1 <= M_in;
	else 
	   begin
            if (!state)
                d1 <= d1 - 1;
            else
                d1 <= d1 + 1;
	   end
end

//======== the counter and S/A block operation ========
always @(posedge gclk or posedge rst)
begin
	if (rst)
			d2 <= 0;
	else
		begin
			// check the counter condition to configure it as up or down counter
			if (!state)				
				d2 <= d2 + d1;				//S/A as adder
			else
				d2 <= d2 - d1 -1;			//S/A as subtractor, extra -1 to generate correct weight
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

assign done = ((d1 == M_in) && (state == 1));

endmodule