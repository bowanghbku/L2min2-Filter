//CIC3 filter runing for 1024 cycles --> max decimation 341   
`timescale 1ns/1ns

module CIC3 (input wire               clk,
			 input wire               rst_in,
			 input wire        [9:0]  M_in,
			 input wire         	  d_in,
			 output wire       [26:0] d_out,
			 output reg               done
			);

// Integrator stage registers
reg [9:0]  d1;
reg [19:0] d2;
reg [26:0] d3;
reg [26:0] d_tmp, d_d_tmp;

// Comb stage registers
reg [26:0] d4, d4_d;
reg [26:0] d5, d5_d;
reg [26:0] d6;
reg [9:0]  count;
reg 	   v_comb;  // Valid signal for comb section running at output rate

//======== Synchronize reset for AD ==========
wire rstb_raw    = !rst_in;
reg  rstb_syn1, rstb;
wire sclk_sync  = clk  && !rstb;
always @(negedge sclk_sync or negedge rstb_raw)
begin
    if (!rstb_raw) begin
        rstb_syn1 <= 1'b0;
        rstb      <= 1'b0;
    end
    else begin
        rstb_syn1 <= 1'b1;
        rstb      <= rstb_syn1;
    end
end    
//======== Synchronize reset for AD ==========
reg [1:0] deci_cnt;
wire tmp_done  = (deci_cnt==3);		//stop the integrator stage after 3 times decimation for sinc3
wire gclk_inte = clk && !tmp_done && rstb;
always @(negedge gclk_inte or negedge rstb)
begin
	if (!rstb)
		begin
			d1 <= 0;
			d2 <= 0;
			d3 <= 0;				
			d_tmp <= 0;
			v_comb <= 0;
			count  <= 0;
			deci_cnt <= 0;
		end 
	else
		begin
			// Integrator section
			d1 <= d_in + d1;				
			d2 <= d1   + d2;				
			d3 <= d3   + d2;
			// Decimation control
			if (count == M_in)
				begin
					d_tmp  	 <= d3;
					v_comb 	 <= 1'b1;
					deci_cnt <= deci_cnt + 1;
					count    <= 1;
				end
			else
				begin
					count <= count + 1;
					v_comb <= 1'b0;
				end
		end
end

wire gclk_comb = clk && !done && rstb && v_comb;
reg [2:0] cnt_comb; //comb stage works for 5 cycles (3 decimation + 2 extra delay) to finish one conversion
always @(negedge gclk_comb or negedge rstb)  //Comb section running at output rate
begin
	if (!rstb)
		begin
			d_d_tmp  <= 0;						
			d4 		 <= 0;
			d4_d 	 <= 0;				
			d5 		 <= 0;
			d5_d 	 <= 0;	
			d6 		 <= 0;	
			cnt_comb <= 0;
			done     <= 0;					
		end
	else
		begin
				// Comb section
				d_d_tmp  <= d_tmp;
				d4 	 	 <= d_tmp 	- d_d_tmp;
				d4_d 	 <= d4;
				d5 	     <= d4 		- d4_d;
				d5_d 	 <= d5;
				d6 	     <= d5 		- d5_d;
				if (cnt_comb == 4)
					done <= 1;
				else
					begin
						cnt_comb <= cnt_comb + 1;
					end
		end
end		
assign d_out = d6;					
endmodule
