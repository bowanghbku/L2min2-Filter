//CIC2 filter runing for 1024 cycles --> max decimation 512   
`timescale 1ns/1ns

module CIC2 (input wire               clk,
			 input wire               rst_in,
			 input wire        [9:0]  M_in,
			 input wire         	  	  d_in,
			 output wire       [18:0] d_out,
			 output reg               done
			);

// Integrator stage registers
reg [9:0]  d1;
reg [18:0] d2;
reg [18:0] d_tmp, d_d_tmp;

// Comb stage registers
reg [18:0] d3, d3_d;
reg [18:0] d4;
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
wire tmp_done  = (deci_cnt==2);		//stop the integrator stage after 2 times decimation for sinc2
wire gclk_inte = clk && !tmp_done && rstb;
always @(negedge gclk_inte or negedge rstb)
begin
	if (!rstb)
		begin
			d1 <= 0;
			d2 <= 0;			
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
			// Decimation control
			if (count == M_in)
				begin
					d_tmp  	 <= d2;
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
reg [1:0] cnt_comb; //comb stage works for 3 cycles (2 decimation + 1 extra delay) to finish one conversion
always @(negedge gclk_comb or negedge rstb)  //Comb section running at output rate
begin
	if (!rstb)
		begin
			d_d_tmp  <= 0;						
			d3 		 <= 0;
			d3_d 	 <= 0;					
			d4 		 <= 0;	
			cnt_comb <= 0;
			done     <= 0;					
		end 
	else
		begin
				// Comb section
				d_d_tmp  <= d_tmp;
				d3 	 	 <= d_tmp 	- d_d_tmp;
				d3_d 	 <= d3;
				d4 	     <= d3 		- d3_d;
				if (cnt_comb == 2)
					done <= 1;
				else
					begin
						cnt_comb <= cnt_comb + 1;
					end
		end
end		
assign d_out = d4;					
endmodule