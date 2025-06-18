import "./Dropdown.css";
import React from "react";
import DropdownButton from "../DropdownButton/DropdownButton";
import DropdownContent from "../DropdownContent/DropdownContent";
import {useState} from "react";
const Dropdown = ({buttonText, contentText}) => {

    const [open, setOpen] = useState(false);
    const toggleDropdown = () => {
        setOpen((open) => !open);
    };

    return <div>
        <DropdownButton toggle={toggleDropdown} open={open}>
            {buttonText}
        </DropdownButton>
        <DropdownContent open={open}>
            {contentText}
        </DropdownContent>
    </div>;
};

export default Dropdown;
